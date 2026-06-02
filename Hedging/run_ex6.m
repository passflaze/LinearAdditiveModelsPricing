function LA_results_es6 = run_ex6(params, market, diagnostics, do_backtest, vega_mode)
%RUN_EX6  Risk Management (Ex.6): delta-gamma-vega hedge of the AB exotics
%         (CoC, PoP, Chooser) with vanilla call, vanilla put and the future.
%
%   Two vega definitions are supported (see the long discussion in
%   greeks_exotic_AB / shock_recalibrate_AB):
%     'scale'   - analytic bump of the normalized scale (self-similar /
%                 sticky-moneyness smile move); (k,eta) held fixed, no
%                 re-calibration, low noise.
%     'recalib' - sticky-strike parallel shift of the implied-vol surface +
%                 AB re-calibration (model-consistent vega); (k,eta) move.
%   Default 'both': computes the two, reports them side by side together with
%   the (k,eta) drift, and BUILDS THE HEDGE with the recalib (sticky-strike)
%   vega when available.
%
% INPUTS:
%   params      - (struct) calibrated parameters; uses params.AB = [k; eta]
%   market      - (struct) curve + surface from run_ex2
%   diagnostics - (optional, default false) verbose prints
%   do_backtest - (optional, default false) recalibrating backtest over the
%                 next two Tuesdays (uses the scale-bump vega for speed)
%   vega_mode   - (optional, default 'both') 'scale' | 'recalib' | 'both'
%
% OUTPUT:
%   LA_results_es6 - (struct) .Greeks .VegaCompare .Positions .Backtest

    addpath("Utilities/");
    addpath("Distributions/");
    addpath(genpath('Pricing/'));
    addpath(genpath('Simulation/'));
    addpath(genpath('Calibration/'));
    addpath(genpath('Hedging/'));

    if nargin < 3 || isempty(diagnostics), diagnostics = false; end
    if nargin < 4 || isempty(do_backtest), do_backtest = false; end
    if nargin < 5 || isempty(vega_mode),   vega_mode   = 'both';  end
    want_recalib = any(strcmpi(vega_mode, {'recalib','both'}));

    fprintf('=========================================================================\n');
    fprintf('        EX.6 - RISK MANAGEMENT: AB EXOTICS DELTA-GAMMA-VEGA HEDGE         \n');
    fprintf('=========================================================================\n\n');

    %% --- AB state on the same (iT1,iT2) window as run_ex4 ----------------
    iT1 = 2;  iT2 = 4;
    I0_AB_val    = I0_AB(0, params.AB);
    sigma_t_AB   = market.sigma_ATM / I0_AB_val;
    scale_factor = [sigma_t_AB(iT1)*sqrt(market.yf(iT1)), ...
                    sigma_t_AB(iT2)*sqrt(market.yf(iT2))];

    F_t0_t2 = market.forward(iT2);

    %% --- Hedge vanilla strikes (Kc != Kp, different moneyness) -----------
    Kc = snap_to_grid(market.strikes, F_t0_t2 + 3);   % near-OTM call
    Kp = snap_to_grid(market.strikes, F_t0_t2 - 7);   % further-OTM put

    mkt = struct('forward', F_t0_t2, ...
                 'K1', 1, ...                 % compound strike (as run_ex4)
                 'K2', F_t0_t2, ...           % inner strike = F(T2,T2) (ATM)
                 'Kc', Kc, 'Kp', Kp, ...
                 'df', [market.discount_factor(iT1), market.discount_factor(iT2)]);

    mc       = struct('N_sim',1e6,'M',16,'dz',5e-3,'N_grid',300,'seed',1234);
    bumps    = struct('dF', 0.5, 'dSig', 1e-2, 'dvol', 0.5);  % dvol: sticky-strike, $/sqrt(yr)
    costRule = struct('fut_bp', 1, 'opt_bp', 4);

    exotics = {'CoC','PoP','Chooser'};
    nE = numel(exotics);

    if diagnostics
        fprintf('  F(t0,T2) = %.4f | hedge strikes  Kc = %.2f  Kp = %.2f\n', F_t0_t2, Kc, Kp);
        fprintf('  scale = [%.5f, %.5f] | bumps dF=%.2f dSig=%.0e dvol=%.2f\n', ...
            scale_factor(1), scale_factor(2), bumps.dF, bumps.dSig, bumps.dvol);
        fprintf('  vega_mode = %s\n\n', vega_mode);
    end

    %% --- Greeks of the hedge instruments (scale-bump, shared) ------------
    gC = greeks_vanilla_AB('call',   params.AB, scale_factor(2), mkt, bumps);
    gP = greeks_vanilla_AB('put',    params.AB, scale_factor(2), mkt, bumps);
    gF = greeks_vanilla_AB('future', params.AB, scale_factor(2), mkt, bumps);

    %% --- Scale-bump exotic greeks ---------------------------------------
    G = repmat(struct('price',[],'delta',[],'gamma',[],'vega',[]), nE, 1);
    for ee = 1:nE
        G(ee) = greeks_exotic_AB(exotics{ee}, params.AB, scale_factor, mkt, mc, bumps);
    end

    %% --- Sticky-strike recalibration vega (shared +/- shocks) -----------
    rec = [];
    if want_recalib
        fprintf('  > Sticky-strike recalibration vega (dvol = %.2f, 2 calibrations)...\n', bumps.dvol);
        rec = vega_recalib_sticky(exotics, market, params.AB, mkt, mc, iT1, iT2, bumps.dvol);
        fprintf('    (k,eta) drift  +dvol: [%+.4f, %+.4f]  -dvol: [%+.4f, %+.4f]\n\n', ...
            rec.drift(1,1), rec.drift(2,1), rec.drift(1,2), rec.drift(2,2));
    end

    %% --- Build the hedge (primary vega = recalib if available) ----------
    P  = repmat(struct('nC',[],'nP',[],'nF',[]), nE, 1);
    initCost = zeros(nE,1);
    resTab   = zeros(nE,3);
    vega_scale_v   = zeros(nE,1);
    vega_recalib_v = nan(nE,1);

    for ee = 1:nE
        gEx = G(ee);
        gCh = gC;  gPh = gP;                 % hedge-instrument greeks
        vega_scale_v(ee) = gEx.vega;

        if want_recalib
            vega_recalib_v(ee) = rec.exotic.(exotics{ee});
            gEx.vega = rec.exotic.(exotics{ee});   % override exotic vega
            gCh.vega = rec.call;                   % consistent vanilla vegas
            gPh.vega = rec.put;
        end

        [P(ee), resG] = build_hedge_AB(gEx, gCh, gPh, gF);
        initCost(ee)  = hedging_cost(P(ee), mkt.forward, costRule);
        resTab(ee,:)  = [resG.delta, resG.gamma, resG.vega];
        if diagnostics
            fprintf('  %-8s | V=%.4f d=%+.4f g=%+.4e | nC=%+.3f nP=%+.3f nF=%+.3f\n', ...
                exotics{ee}, gEx.price, gEx.delta, gEx.gamma, P(ee).nC, P(ee).nP, P(ee).nF);
        end
    end

    %% --- Tables ----------------------------------------------------------
    GreeksT = table(exotics(:), [G.price]', [G.delta]', [G.gamma]', ...
        'VariableNames', {'Exotic','Price','Delta','Gamma'});

    VegaT = table(exotics(:), vega_scale_v, vega_recalib_v, ...
        'VariableNames', {'Exotic','Vega_scale_rel','Vega_recalib_sticky'});

    PosT = table(exotics(:), [P.nC]', [P.nP]', [P.nF]', initCost, ...
        resTab(:,1), resTab(:,2), resTab(:,3), ...
        'VariableNames', {'Exotic','nCall','nPut','nFuture','InitHedgeCost', ...
                          'ResDelta','ResGamma','ResVega'});

    disp(' '); disp('==================== GREEKS: PRICE / DELTA / GAMMA ===================='); disp(GreeksT);
    disp(' '); disp('==================== VEGA: scale-bump vs sticky-strike ==============='); disp(VegaT);
    if want_recalib
        fprintf('  Note: Vega_scale is per RELATIVE vol shift (dSig); Vega_recalib is per\n');
        fprintf('  ABSOLUTE Bachelier-vol unit (dvol). Compare magnitudes only after\n');
        fprintf('  matching units. Hedge built with the sticky-strike (recalib) vega.\n');
    end
    disp(' '); disp('================ HEDGE POSITIONS, COST & RESIDUALS ==============='); disp(PosT);

    %% --- Optional backtest (scale-bump vega for speed) ------------------
    Backtest = struct();
    if do_backtest
        valueDate = market.valueDate;
        tuesdays  = next_tuesdays(valueDate, 2);
        dataOpts  = struct('callpath',"Data/datacalls", 'putpath',"Data/dataputs", ...
                           'expiryFile',"Data/Expiries_Futures.txt");
        for ee = 1:nE
            contract = struct('K1', mkt.K1, 'K2', mkt.K2, 'Kc', Kc, 'Kp', Kp, ...
                              'E1', market.expiries(iT1), 'E2', market.expiries(iT2));
            state0   = struct('params_AB', params.AB, 'scale_factor', scale_factor, 'mkt', mkt);
            try
                fprintf('\n  > Backtest %s over %s ...\n', exotics{ee}, ...
                    strjoin(string(tuesdays,'yyyy-MM-dd'), ', '));
                Bt = hedge_backtest(exotics{ee}, contract, state0, P(ee), ...
                        tuesdays, dataOpts, mc, bumps, costRule);
                disp(Bt);
                Backtest.(exotics{ee}) = Bt;
            catch ME
                fprintf('  [backtest %s skipped] %s\n', exotics{ee}, ME.message);
                Backtest.(exotics{ee}) = ME.message;
            end
        end
    else
        fprintf('\n  (Backtest skipped: call run_ex6(params, market, diag, true) to enable.)\n');
    end

    %% --- Pack ------------------------------------------------------------
    LA_results_es6.Greeks      = GreeksT;
    LA_results_es6.VegaCompare = VegaT;
    LA_results_es6.Positions   = PosT;
    LA_results_es6.Backtest    = Backtest;
    if want_recalib, LA_results_es6.RecalibDiag = rec; end
end

% =========================================================================
function K = snap_to_grid(strikes, target)
    strikes = strikes(:);
    [~, i]  = min(abs(strikes - target));
    K = strikes(i);
end

% =========================================================================
function tue = next_tuesdays(fromDate, n)
    tue = NaT(n,1);
    d = dateshift(fromDate, 'start', 'day');
    cnt = 0;
    while cnt < n
        d = d + caldays(1);
        if weekday(d) == 3   % 1=Sun, 3=Tue
            cnt = cnt + 1;
            tue(cnt) = d;
        end
    end
end
