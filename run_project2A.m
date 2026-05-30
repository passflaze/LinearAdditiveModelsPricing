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
%  I parametri calibrati passano da ex2 a ex3 tramite le struct `params` e
%  `market`.
%  =========================================================================

clear; clc; close all;

% =========================================================================
% PATHS INITIALIZATION
% =========================================================================
addpath("Utilities/");
addpath("ex2/");
addpath("Calibration/");
addpath("Distributions/");
addpath("Functions/");
addpath("Simulation/");
addpath("Simulation/Simulation_MA/");

%% =========================================================================
%  EX 0 - PDF COMPARISON (MA / AB / GL) + AB MONTE CARLO VALIDATION
%  =========================================================================
%  ex0 e uno script di analisi delle distribuzioni: confronta le PDF dei tre
%  modelli additivi e valida la AB via simulazione Monte Carlo.
ex0;

%% =========================================================================
%  EX 1 + EX 2 - CURVE BOOTSTRAP & VOLATILITY SURFACE CALIBRATION
%  =========================================================================
%  calibrate_surface svolge prima l'Esercizio 1 (bootstrap della curva
%  Forward e dei Discount Factor) e poi l'Esercizio 2 (calibrazione della
%  superficie di volatilita). Restituisce:
%    params : parametri calibrati  (params.AB / params.MA / params.GL)
%    market : dati di mercato e supporto (forward, discount_factor,
%             sigma_ATM, yf, expiries, ...)
opts            = struct();
opts.callpath   = "Data/datacalls";
opts.putpath    = "Data/dataputs";
opts.expiryFile = "Data/Expiries_Futures.txt";
opts.valueDate  = datetime(2020, 06, 02);
opts.verbose    = true;    % report completo di calibrazione a video
opts.plot       = true;    % plot delle distribuzioni implicite

[params, market] = calibrate_surface(opts);

%% =========================================================================
%  EX 3 - FORWARD-START PRICING (con i parametri calibrati in ex2)
%  =========================================================================
LA_results = run_ex3(params, market);

%% =========================================================================
%  DONE
%  =========================================================================
fprintf('\n=========================================================================\n');
fprintf('  PROJECT 2A PIPELINE COMPLETED (EX2 calibration -> EX3 forward-start).  \n');
fprintf('=========================================================================\n');
