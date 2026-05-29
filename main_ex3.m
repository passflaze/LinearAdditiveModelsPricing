%% MAIN_EX3  Example caller for the Exercise 3 simulation/pricing engine.
%  The actual "run completo" lives in Simulation/run_ex3_simulation.m, a FUNCTION
%  that takes the already-CALIBRATED model shape parameters as input. This script
%  only collects those parameters and invokes it.
%
%  Calibrated shape parameters (from Exercise 2):
%    - MA: alpha/beta from the MA volatility-surface calibration (objective_function_MA)
%    - GL: alpha/beta from the GL volatility-surface calibration (objective_function_GL)
%    - AB: k/eta from calibrateAB (pooled OTM $-price L2 fit, paper Eq. 20)

clear; clc; close all;

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'Simulation'));
addpath(fullfile(root, 'Utilities'));

% ------------------------------------------------------------------------
% MA and GL: calibrated shape parameters (Exercise 2 outputs)
% ------------------------------------------------------------------------
params.MA = struct('alpha', 1, 'beta', 0.982156);   % Minimal Additive
params.GL = struct('alpha', 0.440992, 'beta', 0.402175);   % Generalized Logistic

% ------------------------------------------------------------------------
% AB: either calibrate now (recommended) or paste your stored Exercise 2 values
% ------------------------------------------------------------------------

params.AB = struct('k', 0.932660, 'eta', -0.063103);

% ------------------------------------------------------------------------
% Run Exercise 3
% ------------------------------------------------------------------------
opts = struct('N_sim', 1e6, 'seed', 2, 'n_strikes', 20, 'doPlots', true);
results = run_ex3_simulation(params, opts);
