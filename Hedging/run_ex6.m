function LA_results_es6 = run_ex6(maturity_index, Kcall, Kput, bump_sigma, CoC_euro, PoP_euro, Ch_euro, products, greeks, tuesdays, dynamic)
% RUN_EX6_HEDGING  Dynamic wrapper for Risk Management Delta-Gamma-Vega hedge.
%
% INPUTS:
%   maturity_index - (struct) Maturities for call, put, future
%   Kcall / Kput   - (scalar or 'ATM') Strikes for the options
%   bump_sigma     - (array) Bumps for finite differences
%   CoC_euro       - (scalar) Budget for Cash-or-Nothing (>0 long, <0 short, 0 exclude)
%   PoP_euro       - (scalar) Budget for Pay-Later     (>0 long, <0 short, 0 exclude)
%   Ch_euro        - (scalar) Budget for Chooser       (>0 long, <0 short, 0 exclude)
%   products       - (cell) Instruments to use (e.g. {'Call', 'Future'})
%   greeks         - (cell) Greeks to hedge (e.g. {'Delta', 'Vega'})
%   tuesdays       - (datetime array) Forward dates for the backtest
%   dynamic        - (logical) True to enable dynamic rebalancing
%
% OUTPUT:
%   LA_results_es6 - (struct) .Greeks .Positions .Cost .Backtest tables

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

    % Disable figure pop-ups for the entire execution to speed it up
    set(0, 'DefaultFigureVisible', 'off');

    fprintf('=========================================================================\n');
    fprintf('        EX.6 - RISK MANAGEMENT: AB EXOTICS DELTA-GAMMA-VEGA HEDGE        \n');
    fprintf('=========================================================================\n\n');

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
    [~, params, market] = evalc('run_ex2(opts)');

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
    %  2. DYNAMIC INPUT RESOLUTION
    %  ========================================================================
    % Parse 'ATM' or empty strikes dynamically
    if isempty(Kcall) || (ischar(Kcall) && strcmpi(Kcall, 'ATM')), Kcall = F_t0_t2; end
    if isempty(Kput)  || (ischar(Kput)  && strcmpi(Kput, 'ATM')),  Kput = F_t0_t2;  end

    params_hedge = struct();
    params_hedge.forward = F_t0_t2;
    params_hedge.K1      = 1;         
    params_hedge.K2      = F_t0_t2;   
    params_hedge.Kc      = Kcall;   
    params_hedge.Kp      = Kput;   

    % Numerical & Simulation Settings
    mc = struct('N_sim', 1e6, 'M', 16, 'dz', 5e-3, 'N_grid', 300, 'seed', 1234);
    % Forward bump for Delta/Gamma: 1e-4 $ was far too small (gamma = second
    % difference / dF^2 ~ /1e-8 -> dominated by FFT-interpolation / MC noise).
    % A ~0.5 $ bump on a ~35-40 $ WTI forward gives stable central-difference
    % Greeks while staying local.
    bumps      = struct('dF', 0.5, 'dSig', bump_sigma);
    costRule   = struct('fut_bp', 1, 'opt_bp', 4);

    %% ========================================================================
    %  3. HEDGING BASKET DEFINITION & INITIAL GREEKS (t0)
    %  ========================================================================
    % A single ATM call+put shares gamma AND vega (put-call parity) and cannot
    % span Delta-Gamma-Vega. We therefore hedge with 3 DISTINCT vanilla strikes
    % on the exotic's terminal leg T2 (ATM + OTM wings) so the 3x3 Greek matrix
    % is full rank. Wing offset = wing_mult * (1 stdev in $).
    %
    % Keep wing_mult MODERATE (~0.5): wings too deep OTM have tiny gamma/vega,
    % which makes the Greek matrix ill-conditioned and blows up the hedge
    % notionals (and the rebalancing cash flows). build_hedge reports cond(A).
    % NB: 'products'/Kput are superseded by this basket (vanillas-only hedge).
    fprintf('Building vanilla hedging basket (3 strikes @ T2)...\n');

    empty_greek = struct('price', 0, 'delta', 0, 'gamma', 0, 'vega', 0);

    vanilla_mat = iT2;
    wing_mult   = 0.5;   % wing distance in stdev units (tune for conditioning)
    dK          = wing_mult * scale_factor_vanilla(vanilla_mat);
    Kc_center   = Kcall;   % resolved ATM anchor (= F if 'ATM' was passed)

    basket = struct( ...
        'kind', {'call',      'call',           'put'}, ...
        'K',    {Kc_center,   Kc_center + dK,   Kc_center - dK}, ...
        'mat',  {vanilla_mat, vanilla_mat,      vanilla_mat});

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
        [~, G(ee)] = evalc('greeks_exotic_AB(active_exotics{ee}, params.AB, scale_factor_exotic, market, params_hedge, mc, bumps)');

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

    % Destroy any hidden plots created in memory
    close all hidden;

    % Re-enable normal plot visibility
    set(0, 'DefaultFigureVisible', 'on');

    fprintf('\nExecution completed successfully. Results packed in LA_results_es6.\n');
end