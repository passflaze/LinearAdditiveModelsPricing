%% =========================================================================
%  MAIN SCRIPT - PROJECT 2A: LINEAR ADDITIVE MODELS (MA / AB / GL)
%  =========================================================================
%  Questo e l'UNICO main del progetto. Orchestra l'intera pipeline:
%
%    EX 0  ->  ex0               : confronto delle PDF dei 3 modelli
%                                  (Asymmetric Laplace MA / AB-NIG / GL) e
%                                  validazione Monte Carlo della AB.
%
%    EX 1  ->  calibrate_surface : bootstrap della curva Forward e dei
%                                  Discount Factor (Esercizio 1 del Project 2),
%                                  svolto all'interno di ex2.
%
%    EX 2  ->  calibrate_surface : bootstrap curva (EX 1) + ATM vol + filtro
%                                  superficie + calibrazione AB / MA / GL.
%                                  Restituisce i PARAMETRI calibrati.
%
%    EX 3  ->  run_ex3           : pricing di opzioni forward-start con i 3
%                                  modelli, usando i parametri calibrati in
%                                  ex2 (nessuna ri-calibrazione, nessun
%                                  parametro hardcoded).
%
%  I parametri calibrati passano da ex2 a ex3 tramite i vettori `params` e
%  `market`.
%  =========================================================================

clear; clc; close all;

% =========================================================================
% PATHS INITIALIZATION
% =========================================================================
addpath("Utilities/");
addpath("Distributions/");
addpath("Calibration/");
addpath("Calibration/Calibration_AB/");
addpath("Calibration/Calibration_MA/");
addpath("Calibration/Calibration_GL/");
addpath("Simulation/");
addpath("Simulation/Simulation_MA/");
addpath("Hedging/");

%% =========================================================================
%  EX 0 - PDF COMPARISON (MA / AB / GL) + AB MONTE CARLO VALIDATION
%  =========================================================================
%  ex0 e uno script di analisi delle distribuzioni: confronta le PDF dei tre
%  modelli additivi e valida la AB via simulazione Monte Carlo.
plot_ex0 = true;
run_ex0(plot_ex0);

%% =========================================================================
%  EX 1 + EX 2 - CURVE BOOTSTRAP & VOLATILITY SURFACE CALIBRATION
%  =========================================================================
%  calibrate_surface svolge prima l'Esercizio 1 (bootstrap della curva
%  Forward e dei Discount Factor) e poi l'Esercizio 2 (calibrazione della
%  superficie di volatilita). Restituisce:
%    params : parametri calibrati come vettori colonna
%             (params.AB = [k;eta], params.MA = [alpha;beta], params.GL = [alpha;beta])
%    market : dati di mercato e supporto (forward, discount_factor,
%             sigma_ATM, yf, expiries, ...)
opts            = struct();
opts.callpath   = "Data/datacalls";
opts.putpath    = "Data/dataputs";
opts.expiryFile = "Data/Expiries_Futures.txt";
opts.valueDate  = datetime(2020, 06, 02);
opts.verbose    = false;    % report completo di calibrazione a video
opts.plot       = false;    % plot delle distribuzioni implicite

[params, market] = run_ex2(opts);   %%DA CAMBIARE METTIAMO SOLO IL PRINT DELLA CHIAMATA ALLA FUNZIONE E DEI PARAMETRI FINALI 


%% =========================================================================
%  EX 3 - FORWARD-START PRICING (con i parametri calibrati in ex2)
%  =========================================================================
opts_ex3 = struct();
opts_ex3.plot = false;
opts_ex3.verbose = false;

LA_results = run_ex3(params, market, opts_ex3);

%% =========================================================================
%  EX 4 - CoC-PoP-Chooser PRICING (con i parametri calibrati in ex2)
%  =========================================================================
opts_ex4 = struct();
opts_ex4.plot = false;
opts_ex4.verbose = false;
opts_ex4.smart_extrap = false;

LA_results_es4 = run_ex4(params, market, opts_ex4);


%% =========================================================================
%  EX 6 - HEDGIING 
%  =========================================================================
clc;

maturity_index = struct();
maturity_index.call = 4;
maturity_index.put = 4;
maturity_index.future = 4;

Kcall = 'ATM';
Kput = 'ATM'; % check if ATM is dynamic
% Vol bump for the (recalibrated) vega: 1e-4 = 1 bp was below the fmincon
% convergence tolerance -> vega dominated by optimizer noise. 1e-2 = 1 vol
% point gives a clean central-difference signal.
bump_sigma = 1e-2 * ones(8,1);

CoC_euro = -1e6; 
PoP_euro = 1e4; 
Ch_euro  = -1e6;

% The hedge basket is now 3 distinct vanilla strikes @ T2 (built inside
% run_ex6), which spans Delta-Gamma-Vega. 'products' is kept for signature
% compatibility but superseded by the vanilla basket.
products = {'Call', 'Put'};
greeks   = {'Delta', 'Gamma', 'Vega'};
tuesdays = [datetime(2020, 6, 9); datetime(2020, 6, 16)];
dynamic = true;
LA_results_es6 = run_ex6(maturity_index, Kcall, Kput, bump_sigma, CoC_euro, PoP_euro, Ch_euro, products, greeks, tuesdays, dynamic);

%% =========================================================================
%  DONE
%  =========================================================================
fprintf('\n=========================================================================\n');
fprintf('  PROJECT 2A PIPELINE COMPLETED (EX2 calibration -> EX3 forward-start -> EX4 CoC-PoP-Chooser pricing).  \n');
fprintf('=========================================================================\n');
