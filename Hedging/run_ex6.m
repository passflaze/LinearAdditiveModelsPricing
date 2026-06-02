function LA_results_es6 = run_ex6(params, market, diagnostics)
%RUN_EX6  Risk Management (Ex.6): delta-gamma-vega hedge of the AB exotics
%         (CoC, PoP, Chooser) with vanilla call, vanilla put and the front
%         future, plus a P&L / cost backtest over the next two Tuesdays.
%
%   Mirrors the structure of run_ex4: consumes the params/market produced by
%   ex2 (no re-calibration, no hardcoded parameters) and orchestrates the
%   Hedging/ module.
%
% INPUTS:
%   params      - (struct) calibrated parameters; uses params.AB = [k; eta]
%   market      - (struct) curve + surface (forward, discount_factor,
%                          sigma_ATM, yf, expiries, strikes, ...)
%   diagnostics - (optional, default false) verbose prints
%
% OUTPUT:
%   LA_results_es6 - (struct) .Greeks .Positions .Cost .Backtest tables

    addpath("Utilities/");
    addpath("Distributions/");
    addpath(genpath('Pricing/'));
    addpath(genpath('Simulation/'));
    addpath(genpath('Calibration/'));
    addpath(genpath('Hedging/'));

    if nargin < 3 || isempty(diagnostics), diagnostics = false; end

    fprintf('=========================================================================\n');
    fprintf('        EX.6 - RISK MANAGEMENT: AB EXOTICS DELTA-GAMMA-VEGA HEDGE         \n');
    fprintf('=========================================================================\n\n');

    %% --- Unpack market & AB state (same window as run_ex4: iT1=2, iT2=4) --
    iT1 = 2;  iT2 = 4;
    I0_AB_val    = I0_AB(0, params.AB);
    sigma_t_AB   = market.sigma_ATM / I0_AB_val;
    scale_factor = [sigma_t_AB(iT1)*sqrt(market.yf(iT1)), ...
                    sigma_t_AB(iT2)*sqrt(market.yf(iT2))];

    F_t0_t2 = market.forward(iT2);
    mkt = struct();
    mkt.forward = F_t0_t2;
    mkt.K1      = 1;                 % compound strike (as run_ex4)
    mkt.K2      = F_t0_t2;           % inner strike = F(T2,T2)  (ATM)
    mkt.df      = [market.discount_factor(iT1), market.discount_factor(iT2)];
    mkt.Kc      = F_t0_t2;           % hedge call strike  (TODO: choose grid strike)
    mkt.Kp      = F_t0_t2;           % hedge put  strike  (TODO: choose grid strike)

    mc       = struct('N_sim',1e6,'M',16,'dz',5e-3,'N_grid',300,'seed',1234);
    bumps    = struct('dF', 0.5, 'dSig', 1e-2);
    costRule = struct('fut_bp', 1, 'opt_bp', 4);

    exotics = {'CoC','PoP','Chooser'};

    %% --- Greeks + hedge per exotic ---------------------------------------
    gC = greeks_vanilla_AB('call',   params.AB, scale_factor(2), mkt, bumps);
    gP = greeks_vanilla_AB('put',    params.AB, scale_factor(2), mkt, bumps);
    gF = greeks_vanilla_AB('future', params.AB, scale_factor(2), mkt, bumps);

    nE = numel(exotics);
    G  = repmat(struct('price',[],'delta',[],'gamma',[],'vega',[]), nE, 1);
    P  = repmat(struct('nC',[],'nP',[],'nF',[]), nE, 1);
    initCost = zeros(nE,1);

    for ee = 1:nE
        G(ee) = greeks_exotic_AB(exotics{ee}, params.AB, scale_factor, mkt, mc, bumps);
        [P(ee), resG] = build_hedge_AB(G(ee), gC, gP, gF);
        initCost(ee)  = hedging_cost(P(ee), mkt.forward, costRule);
        if diagnostics
            fprintf('  %-8s | dV=%+.4f dG=%+.4e nC=%+.3f nP=%+.3f nF=%+.3f | resid d/g/v = %.1e/%.1e/%.1e\n', ...
                exotics{ee}, G(ee).delta, G(ee).gamma, P(ee).nC, P(ee).nP, P(ee).nF, ...
                resG.delta, resG.gamma, resG.vega);
        end
    end

    %% --- Tables ----------------------------------------------------------
    GreeksT = table(exotics(:), [G.price]', [G.delta]', [G.gamma]', [G.vega]', ...
        'VariableNames', {'Exotic','Price','Delta','Gamma','Vega'});
    PosT = table(exotics(:), [P.nC]', [P.nP]', [P.nF]', initCost, ...
        'VariableNames', {'Exotic','nCall','nPut','nFuture','InitHedgeCost'});

    disp(' '); disp('===================== GREEKS (AB) ====================='); disp(GreeksT);
    disp(' '); disp('================ HEDGE POSITIONS & COST ==============='); disp(PosT);

    %% --- Backtest over the next two Tuesdays ------------------------------
    % TODO(Persona B): set the real value date and Tuesdays + data paths, then
    % enable. Requires recalibrate_AB inside hedge_backtest to be wired.
    tuesdays = [datetime(2020,6,9); datetime(2020,6,16)];
    dataOpts = struct('callpath',"Data/datacalls", 'putpath',"Data/dataputs", ...
                      'expiryFile',"Data/Expiries_Futures.txt", 'iT1',iT1, 'iT2',iT2);

    Backtest = struct();
    for ee = 1:nE
        state0 = struct('params_AB', params.AB, 'scale_factor', scale_factor, 'mkt', mkt);
        try
            Backtest.(exotics{ee}) = hedge_backtest(exotics{ee}, state0, P(ee), ...
                tuesdays, dataOpts, mc, bumps, costRule);
        catch ME
            fprintf('  [backtest %s skipped] %s\n', exotics{ee}, ME.message);
            Backtest.(exotics{ee}) = ME.message;
        end
    end

    %% --- Pack ------------------------------------------------------------
    LA_results_es6.Greeks    = GreeksT;
    LA_results_es6.Positions = PosT;
    LA_results_es6.Backtest  = Backtest;
end
