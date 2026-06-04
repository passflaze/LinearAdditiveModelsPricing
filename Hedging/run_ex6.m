function LA_results_es6 = run_ex6(products, hedge_strike, hedge_mat, bump_sigma, CoC_euro, PoP_euro, Ch_euro, greeks, tuesdays, dynamic)
% RUN_EX6_HEDGING  Dynamic wrapper for Risk Management Delta-Gamma-Vega hedge.
%
%   The hedging basket is driven ENTIRELY by the caller: one instrument per
%   entry of 'products'. Examples:
%       {'Call','Put'}           -> 1 call + 1 put
%       {'Call','Call','Put'}    -> 2 calls (same/different K) + 1 put
%       {'Call','Put','Future'}  -> 1 call + 1 put + 1 future
%   Strikes and maturities for each leg are also chosen by the caller via
%   hedge_strike and hedge_mat (see below). For an EXACT hedge use as many
%   independent instruments as targeted Greeks (e.g. 3 for Delta-Gamma-Vega);
%   same-strike call & put are gamma/vega-redundant (build_hedge warns).
%
% INPUTS:
%   products        - (cell) one kind per instrument: 'Call' | 'Put' | 'Future'
%   hedge_strike - (vector) per-leg strike as a $-offset from the ATM forward
%                     at that leg (0 = ATM). Ignored for 'Future'. Same length
%                     as products.
%   hedge_mat       - (vector) per-leg maturity index (e.g. 4 = T2). Same length
%                     as products.
%   bump_sigma      - (array) Bumps for finite differences (vol)
%   CoC_euro        - (scalar) Budget for Cash-or-Nothing (>0 long, <0 short, 0 exclude)
%   PoP_euro        - (scalar) Budget for Pay-Later     (>0 long, <0 short, 0 exclude)
%   Ch_euro         - (scalar) Budget for Chooser       (>0 long, <0 short, 0 exclude)
%   greeks          - (cell) Greeks to hedge (e.g. {'Delta','Gamma','Vega'})
%   tuesdays        - (datetime array) Forward dates for the backtest
%   dynamic         - (logical) True to enable dynamic rebalancing
%
% OUTPUT:
%   LA_results_es6 - (struct) .Greeks .Basket .Weights .Residuals .Cost .Backtest

    %% ========================================================================
    %  0. ENVIRONMENT SETUP
    %  ========================================================================
    % Addpaths usually done in your main folder structure
    addpath("Utilities/");
    addpath("Distributions/");
    addpath(genpath('Pricing/'));
    addpath(genpath('Simulation/'));
    addpath(genpath('Calibration/'));
    addpath(genpath('Hedging/'));


    fprintf('=========================================================================\n');
    fprintf('        EX.6 - RISK MANAGEMENT: AB EXOTICS DELTA-GAMMA-VEGA HEDGE        \n');
    fprintf('=========================================================================\n\n');

    %% ========================================================================
    %  0b. HEDGE FEASIBILITY CHECK (fail fast, before any pricing)
    %  ========================================================================
    % An exact hedge of M Greeks needs at least M independent instruments.
    % With fewer instruments the system is under-determined: build_hedge would
    % silently return a least-squares fit that does NOT zero all the targeted
    % Greeks. Block here with a clear message instead.
    nInstr  = numel(products);
    nGreeks = numel(greeks);
    if nInstr < nGreeks
        warning('run_ex6:UnderHedged', ...
            'Under-determined hedge: %d instruments for %d Greeks.', nInstr, nGreeks);
        error('run_ex6:UnderHedged', ...
            ['Cannot hedge {%s} (%d Greeks) with only %d instrument(s) {%s}. ', ...
             'Add instruments so numel(products) >= numel(greeks) (with DISTINCT ', ...
             'strikes/maturities to keep them independent).'], ...
            strjoin(greeks, ','), nGreeks, nInstr, strjoin(products, ','));
    end

    %% ========================================================================
    %  1. MARKET CALIBRATION & STATE UNPACKING (t0)
    %  ========================================================================
    opts            = struct();
    opts.callpath   = "Data/datacalls";
    opts.putpath    = "Data/dataputs";
    opts.expiryFile = "Data/Expiries_Futures.txt";
    opts.valueDate  = datetime(2020, 06, 02);
    opts.verbose    = false;    
    opts.plot       = false;    

    fprintf('Calibrating market data...\n');
    [params, market] = run_ex2(opts);

    % Time indices
    iT1 = 2;  
    iT2 = 4;

    % Compute scale factors
    I0_AB_val = I0_AB(0, params.AB);
    sigma_t_AB = market.sigma_ATM / I0_AB_val;
    scale_factor_exotic  = [sigma_t_AB(iT1)*sqrt(market.yf(iT1)), ...
                            sigma_t_AB(iT2)*sqrt(market.yf(iT2))];
    scale_factor_vanilla = sigma_t_AB .* sqrt(market.yf);
    F_t0_t2 = market.forward(iT2);

    %% ========================================================================
    %  2. PRICING PARAMETERS
    %  ========================================================================
    params_hedge = struct();
    params_hedge.forward = F_t0_t2;
    params_hedge.K1      = 1;
    params_hedge.K2      = F_t0_t2;
    params_hedge.Kc      = F_t0_t2;   % ATM default (only used by the exotic vega path)

    % Numerical & Simulation Settings
    mc = struct('N_sim', 1e6, 'M', 16, 'dz', 5e-3, 'N_grid', 300, 'seed', 1234);
    % Forward bump for Delta/Gamma: 1e-4 $ was far too small (gamma = second
    % difference / dF^2 ~ /1e-8 -> dominated by FFT-interpolation / MC noise).
    % A ~0.5 $ bump on a ~35-40 $ WTI forward gives stable central-difference
    % Greeks while staying local.
    bumps      = struct('dF', 1e-4*F_t0_t2, 'dSig', bump_sigma);
    costRule   = struct('fut_bp', 1, 'opt_bp', 4);

    %% ========================================================================
    %  3. HEDGING BASKET (driven by the caller's 'products' spec) & GREEKS (t0)
    %  ========================================================================
    % One instrument per entry of 'products'. Strike of each leg is the ATM
    % forward at its maturity plus the caller's $-offset (hedge_strike);
    % 'Future' legs ignore the offset. For an exact Delta-Gamma-Vega hedge use
    % >= 3 independent instruments (distinct strikes); a same-strike call+put
    % share gamma & vega (build_hedge warns via rank / cond(A)).
    fprintf('Building hedging basket from products spec...\n');

    empty_greek = struct('price', 0, 'delta', 0, 'gamma', 0, 'vega', 0);

    nP = numel(products);
    if numel(hedge_strike) ~= nP || numel(hedge_mat) ~= nP
        error('run_ex6:basketSpec', ...
              ['products (%d), hedge_strike (%d) and hedge_mat (%d) must ', ...
               'have the same length.'], nP, numel(hedge_strike), numel(hedge_mat));
    end

    % Target strike = ATM forward(leg) + $ offset, then SNAP to the nearest
    % quoted strike so the hedge instruments are actually tradable. Price and
    % Greeks stay AB (model-consistent with the exotic); at a calibrated strike
    % the AB price ~ market mid.
    basket = repmat(struct('kind', '', 'K', NaN, 'mat', NaN), nP, 1);
    for j = 1:nP
        kind = lower(products{j});
        matj = hedge_mat(j);
        if strcmpi(kind, 'future')
            Kj = NaN;                                   
        else
            if iscell(hedge_strike)
                hs_val = hedge_strike{j};
            else
                hs_val = hedge_strike(j);
            end
            
            if isnumeric(hs_val) && hs_val == 0
                Kj = market.forward(matj);
            else
                Kj = hs_val;
            end
        end
        basket(j) = struct('kind', kind, 'K', Kj, 'mat', matj);
    end

    % Warn if snapping collapsed two vanilla legs onto the SAME (kind,K,mat):
    % they would be identical instruments -> rank-deficient hedge.
    vkey = arrayfun(@(b) sprintf('%s_%.6g_%d', b.kind, b.K, b.mat), basket, 'uni', 0);
    if numel(unique(vkey)) < nP
        warning('run_ex6:DuplicateStrikes', ...
            ['Snapping collapsed two legs onto the same quoted strike: the ', ...
             'basket has redundant instruments. Spread the hedge_moneyness ', ...
             'values further apart.']);
    end

    fprintf('Computing initial Greeks for the hedging basket...\n');
    instr = price_hedge_basket(basket, params.AB, scale_factor_vanilla, market, params_hedge, mc, bumps);

    % Filter active exotics by NON-ZERO budget. Sign convention:
    %   budget > 0  -> LONG  ~ floor(+budget/price) units
    %   budget < 0  -> SHORT ~ -floor(|budget|/price) units
    %   budget = 0  -> excluded
    fprintf('Computing initial Greeks for ACTIVE Exotic Portfolio...\n');
    all_exotics = {'CoC', 'PoP', 'Chooser'};
    all_budgets = [CoC_euro, PoP_euro, Ch_euro];

    active_idx = all_budgets ~= 0;
    active_exotics = all_exotics(active_idx);
    active_budgets = all_budgets(active_idx);
    nE = numel(active_exotics);

    G = repmat(empty_greek, max(1, nE), 1);
    w_exotics = zeros(nE, 1);

    port_t0_price = 0; port_delta = 0; port_gamma = 0; port_vega = 0;

    for ee = 1:nE
        G(ee) = greeks_exotic_AB(active_exotics{ee}, params.AB, scale_factor_exotic, market, params_hedge, mc, bumps);

        % Signed quantity: |budget| sets the size, sign sets long/short.
        w_exotics(ee) = sign(active_budgets(ee)) * floor(abs(active_budgets(ee)) / G(ee).price);

        % Aggregate unitary portfolio Greeks
        port_t0_price = port_t0_price + w_exotics(ee) * G(ee).price;
        port_delta    = port_delta    + w_exotics(ee) * G(ee).delta;
        port_gamma    = port_gamma    + w_exotics(ee) * G(ee).gamma;
        port_vega     = port_vega     + w_exotics(ee) * G(ee).vega;
    end

    %% ========================================================================
    %  4. STATIC HEDGING RESOLUTION & TRANSACTION COSTS (t0)
    %  ========================================================================
    portfolio_greeks.delta = port_delta;
    portfolio_greeks.gamma = port_gamma;
    portfolio_greeks.vega  = port_vega;

    [w0, risk_res] = build_hedge(portfolio_greeks, instr, greeks);

    % Initial transaction costs (slippage + premium) of opening the hedge.
    costs = hedging_cost(w0, instr, costRule);
    fprintf('  Total Initial Hedging Cost: $%+.4e\n\n', costs.total);

    %% ========================================================================
    %  5. DYNAMIC BACKTESTING OVER FUTURE DATES
    %  ========================================================================
    fprintf('Initializing dynamic backtest...\n');
    
    % Build initial state at t0. cost_t0 = day-0 SLIPPAGE only: the premium
    % paid to open the hedge is already embedded in V_hedge (self-financing
    % baseline), so only the bid-ask friction is a true cost.
    prices0 = [instr.price]';
    state_t0 = struct( ...
        'V_exotic', port_t0_price, ...
        'w',        w0, ...
        'V_hedge',  w0' * prices0, ...
        'cost_t0',  costs.slippage ...
    );

    % Pack pricing parameters required inside the forward loop.
    % Non-scalar values (w_exotics, basket, active_exotics) are wrapped in a
    % cell so struct() does NOT replicate prc_params into a struct array.
    prc_params = struct( ...
        'exotics',      {active_exotics}, ...
        'w_exotics',    {w_exotics}, ...
        'params_hedge', params_hedge, ...
        'mc',           mc, ...
        'bumps',        bumps, ...
        'iT1',          iT1, ...
        'iT2',          iT2, ...
        'basket',       {basket} ...
    );

    % Define behavior for dynamic hedging. The rebalancing trigger is RELATIVE:
    % rebalance a Greek only when the residual exceeds rel_tol * |exotic Greek|
    % (i.e. the hedge error exceeds 5% of the gross exposure). Absolute 0.05
    % thresholds vs Greeks of order 1e3-1e5 made the trigger fire every step.
    hedge_rules = struct( ...
        'is_dynamic', dynamic, ...
        'rel_tol',    0.05, ...
        'costRule',   costRule ...
    );

    opts_array = repmat(opts, length(tuesdays), 1);
    for t = 1:length(tuesdays)
        opts_array(t).valueDate = tuesdays(t);
    end

    % Execute the backtest
    try
        Backtest_Results = hedge_backtest(state_t0, opts_array, prc_params, hedge_rules, greeks);
    catch ME
        fprintf('\n[Backtest Failed] %s\n', ME.message);
        Backtest_Results = ME.message;
    end

    %% ========================================================================
    %  6. RESULTS PACKAGING & CLEANUP
    %  ========================================================================
    LA_results_es6 = struct();
    LA_results_es6.Greeks    = G;          % exotic per-type Greeks
    LA_results_es6.Basket    = instr;      % hedge instruments + their Greeks
    LA_results_es6.Weights   = w0;         % static hedge quantities (Nx1)
    LA_results_es6.Residuals = risk_res;   % residual targeted Greeks at t0
    LA_results_es6.Cost      = costs.total;
    LA_results_es6.Backtest  = Backtest_Results;



    fprintf('\nExecution completed successfully. Results packed in LA_results_es6.\n');
end