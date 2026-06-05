%% =========================================================================
%  MAIN SCRIPT - PROJECT 2A: LINEAR ADDITIVE MODELS (MA / AB / GL)
%  =========================================================================
%  This is the ONLY main script of the project. It acts as a centralized
%  orchestrator for the entire pipeline, executing calibration, pricing,
%  sanity checks, and dynamic hedging using three Linear Additive Models:
%  Additive Bachelier (AB), Minimal Additive (MA), and Generalized Logistic (GL).
%
%  PIPELINE OVERVIEW & FUNCTION OUTPUTS:
%
%  EX 0 -> run_ex0:
%     Visual comparison of the Probability Density Functions (PDF) of the
%     3 models (Asymmetric Laplace MA, AB-NIG, GL) and Monte Carlo
%     validation of the AB increment generation. 
%
%  EX 1 & 2 -> run_ex2:
%     Performs the bootstrap of the Forward curve and Discount Factors
%     (Ex. 1), applies a liquidity filter to the volatility surface,
%     and calibrates the AB, MA, and GL parameters (Ex. 2).
%     OUTPUTS:
%       - params : Struct with calibrated parameters (column vectors).
%                  .AB -> [k; eta]
%                  .MA -> [alpha; beta]
%                  .GL -> [alpha; beta]
%       - market : Struct with market data and curves (.strikes, .calls,
%                  .puts, .expiries, .valueDate, .discount_factor, .forward,
%                  .sigma_ATM, .yf, .c_ATM, .moneyness_modified, .R2).
%
%  EX 3 -> run_ex3:
%     Prices forward-start options using the calibrated parameters.
%     OUTPUTS:
%       - LA_results_es3 : Struct containing `.PricingData`. This is a table
%                          comparing Analytic, Monte Carlo, and Uniform Lewis
%                          prices (for ATM), or Monte Carlo vs Confidence
%                          Intervals (for K2 != ATM), plus differences in bps.
%
%  EX 4 -> run_ex4:
%     Prices exotic path-dependent options (Call-on-Call, Put-on-Put, Chooser).
%     OUTPUTS:
%       - LA_results_es4 : Struct containing:
%           .Pricing            -> Table with MC vs Analytic prices, errors, and CIs.
%           .SanityChecks       -> Table validating boundary conditions (e.g., K1=0).
%           .SmartExtrapolation -> Table evaluating CDF truncation effects.
%
%  EX 6 -> run_ex6 & plot_strategy_comparison:
%     Executes dynamic hedging (Delta, Gamma, Vega) on exotic portfolios
%     across multiple time steps (backtesting).
%     OUTPUTS:
%       - LA_results_es6 : Struct containing 6 sub-structs (test1 to test6).
%                          Each contains:
%                          .Greeks    (Exotic per-type Greeks)
%                          .Basket    (Hedge instruments + their Greeks)
%                          .Weights   (Static hedge quantities)
%                          .Residuals (Residual targeted Greeks at t0)
%                          .Cost      (Total initial transaction cost)
%                          .Backtest  (Struct with step-by-step arrays of PnL,
%                                      costs, and residual Greeks over time).
%
%  NOTE: Calibrated parameters flow from ex2 to ex3, ex4, and ex6 seamlessly 
%  via the `params` and `market` structures (no hardcoded parameters).
%  =========================================================================

clear; clc; close all;

addpath("Run_Functions/");
addpath("Utilities/");
addpath("Distributions/");
addpath("Calibration/");
addpath("Calibration/Calibration_AB/");
addpath("Calibration/Calibration_MA/");
addpath("Calibration/Calibration_GL/");
addpath("Simulation/");
addpath("Simulation/Simulation_MA/");
addpath("Hedging/");

%% EX 0 — PDF comparison
plot_ex0 = false;
run_ex0(plot_ex0);

%% EX 1+2 — Curve bootstrap and surface calibration
opts_ex2            = struct();
opts_ex2.callpath   = "Data/datacalls";
opts_ex2.putpath    = "Data/dataputs";
opts_ex2.expiryFile = "Data/Expiries_Futures.txt";
opts_ex2.valueDate  = datetime(2020, 06, 02);

opts_ex2.verbose    = false;   
opts_ex2.plot       = false;    

[params, market] = run_ex2(opts_ex2);   

%% EX 3 — Forward-start pricing
opts_ex3 = struct();
opts_ex3.plot = false;
opts_ex3.verbose = false;
opts_ex3.K2 = 1;

opts_ex3.mc = struct();
opts_ex3.mc.dz_MA = 5e-3;
opts_ex3.mc.Nsim_MA = 0;
opts_ex3.mc.M_MA = 16;
opts_ex3.mc.dz_GL = 5e-3;
opts_ex3.mc.Nsim_GL = 0;
opts_ex3.mc.M_GL = 16;
opts_ex3.mc.dz_AB = 5e-2;
opts_ex3.mc.Nsim_AB = 0;
opts_ex3.mc.M_AB = 16;

LA_results_es3 = run_ex3(params, market, opts_ex3);

%% EX 4 — Exotic pricing (CoC, PoP, Chooser)
opts_ex4 = struct();
opts_ex4.plot = false;
opts_ex4.verbose = false;
opts_ex4.smart_extrap = false;

opts_ex4.K1 = 1;
opts_ex4.K2 = 'ATM';
opts_ex4.mc = struct();
opts_ex4.mc.dz_MA = 5e-3;
opts_ex4.mc.Nsim_MA = 0;
opts_ex4.mc.M_MA = 16;
opts_ex4.mc.dz_GL = 5e-3;
opts_ex4.mc.Nsim_GL = 0;
opts_ex4.mc.M_GL = 16;
opts_ex4.mc.dz_AB = 5e-2;
opts_ex4.mc.Nsim_AB = 0;
opts_ex4.mc.M_AB = 20;

LA_results_es4 = run_ex4(params, market, opts_ex4);

%% EX 6 — Dynamic hedging

%% Scenario 1 — Long Chooser, Delta-Vega hedge
bump_sigma = 1e-2 * ones(8,1);
CoC_euro = 0; 
PoP_euro = 0; 
Ch_euro  = 1e6;

% --- HEDGING BASKET (fully defined here) ---------------------------------
products        = {'Call', 'Put'};
hedge_strike    = [   0,      0];      % 0 means ATM , otherwise specify the strike
hedge_mat       = [   4,      2];
greeks          = {'Delta','Vega'};
tuesdays = [datetime(2020, 6, 9); datetime(2020, 6, 16)];
dynamic  = true;

fprintf('\n--- TEST 1: Long Chooser Delta-Vega Hedging ---\n');
fprintf('    Hedge Basket: 1 Call, 1 Put\n');
fprintf('    Target Greeks: Delta, Vega\n');
LA_results_es6.test1 = run_ex6(products, hedge_strike, hedge_mat, bump_sigma, ...
                         CoC_euro, PoP_euro, Ch_euro, greeks, tuesdays, dynamic, params, market);

%% Scenario 2 — Long Chooser, Gamma-Vega hedge
bump_sigma = 1e-2 * ones(8,1);
CoC_euro = 0; 
PoP_euro = 0; 
Ch_euro  = 1e6;

% --- HEDGING BASKET (fully defined here) ---------------------------------
products        = {'Call', 'Put'};
hedge_strike    = [   0,      0];
hedge_mat       = [   4,      2];
greeks          = {'Gamma','Vega'};
tuesdays = [datetime(2020, 6, 9); datetime(2020, 6, 16)];
dynamic  = true;

fprintf('\n--- TEST 2: Long Chooser Gamma-Vega Hedging ---\n');
fprintf('    Hedge Basket: 1 Call, 1 Put\n');
fprintf('    Target Greeks: Gamma, Vega\n');
LA_results_es6.test2 = run_ex6(products, hedge_strike, hedge_mat, bump_sigma, ...
                         CoC_euro, PoP_euro, Ch_euro, greeks, tuesdays, dynamic, params, market);

%% Scenario 3 — Long CoC+PoP, Gamma-Vega hedge
bump_sigma = 1e-2 * ones(8,1);
CoC_euro = 1e6; 
PoP_euro = 1e6; 
Ch_euro  = 0;

% --- HEDGING BASKET (fully defined here) ---------------------------------
products        = {'Call', 'Put'};
hedge_strike    = [   0,      0];
hedge_mat       = [   4,      2];
greeks          = {'Gamma','Vega'};
tuesdays = [datetime(2020, 6, 9); datetime(2020, 6, 16)];
dynamic  = true;

fprintf('\n--- TEST 3: Long CoC, PoP Gamma-Vega Hedging ---\n');
fprintf('    Hedge Basket: 1 Call, 1 Put\n');
fprintf('    Target Greeks: Gamma, Vega\n');
LA_results_es6.test3 = run_ex6(products, hedge_strike, hedge_mat, bump_sigma, ...
                         CoC_euro, PoP_euro, Ch_euro, greeks, tuesdays, dynamic, params, market);

%% Scenario 4 — Long CoC+PoP, Delta-Gamma-Vega hedge
bump_sigma = 1e-2 * ones(8,1);
CoC_euro = 1e6; 
PoP_euro = 1e6; 
Ch_euro  = 0;

% --- HEDGING BASKET (fully defined here) ---------------------------------
products        = {'Call', 'Put', 'Future'};
hedge_strike    = [   0,      0,         0];
hedge_mat       = [   4,      2,         4];
greeks          = {'Gamma','Vega', 'Delta'};
tuesdays = [datetime(2020, 6, 9); datetime(2020, 6, 16)];
dynamic  = true;

fprintf('\n--- TEST 4: Long CoC, PoP Delta-Gamma-Vega Hedging ---\n');
fprintf('    Hedge Basket: 1 Call, 1 Put, 1 Future\n');
fprintf('    Target Greeks: Gamma, Vega, Delta\n');
LA_results_es6.test4 = run_ex6(products, hedge_strike, hedge_mat, bump_sigma, ...
                         CoC_euro, PoP_euro, Ch_euro, greeks, tuesdays, dynamic, params, market);

%% Scenario 5 — Long CoC + Short Chooser, Vega hedge
bump_sigma = 1e-2 * ones(8,1);
CoC_euro = 1e6; 
PoP_euro = 0; 
Ch_euro  = -1e6;

% --- HEDGING BASKET (fully defined here) ---------------------------------
products        = {'Call'};
hedge_strike    = [   0];
hedge_mat       = [   4];
greeks          = {'Vega'};
tuesdays = [datetime(2020, 6, 9); datetime(2020, 6, 16)];
dynamic  = true;

fprintf('\n--- TEST 5: Long CoC, Short Chooser Vega Hedging ---\n');
fprintf('    Hedge Basket: 1 Call\n');
fprintf('    Target Greeks: Vega\n');
LA_results_es6.test5 = run_ex6(products, hedge_strike, hedge_mat, bump_sigma, ...
                         CoC_euro, PoP_euro, Ch_euro, greeks, tuesdays, dynamic, params, market);

%% Scenario 6 — Long CoC + Short Chooser, Delta hedge
bump_sigma = 1e-2 * ones(8,1);
CoC_euro = 1e6; 
PoP_euro = 0; 
Ch_euro  = -1e6;

% --- HEDGING BASKET (fully defined here) ---------------------------------
products        = {'Call'};
hedge_strike    = [ 0];
hedge_mat       = [ 4];
greeks          = {'Delta'};
tuesdays = [datetime(2020, 6, 9); datetime(2020, 6, 16)];
dynamic  = true;

fprintf('\n--- TEST 6: Long CoC, Short Chooser Delta Hedging ---\n');
fprintf('    Hedge Basket: 1 Call\n');
fprintf('    Target Greeks: Delta\n');
LA_results_es6.test6 = run_ex6(products, hedge_strike, hedge_mat, bump_sigma, ...
                         CoC_euro, PoP_euro, Ch_euro, greeks, tuesdays, dynamic, params, market);


%% -------------------------------------------------------------------------------

fprintf('\n=========================================================================\n');
fprintf('  PROJECT 2A PIPELINE COMPLETED\n');
fprintf('=========================================================================\n');