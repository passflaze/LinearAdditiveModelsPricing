%% MAIN SCRIPT - VOLATILITY SURFACE CALIBRATION (AB, MA & GL MODELS)
% This script loads market option data, bootstraps the forward curve and
% discount factors, calibrates the ATM volatility, and fits both the
% Minimal Additive (MA) and Generalized Laplace (GL) models to the market
% surface.
%
% La logica di calibrazione e ora incapsulata nella funzione riutilizzabile
% calibrate_surface (ex2/calibrate_surface.m), che restituisce i parametri
% calibrati. Questo script e un semplice wrapper che la invoca con i default,
% stampa il report e plotta le distribuzioni implicite.

clear; clc; close all;

% =========================================================================
% PATHS INITIALIZATION
% =========================================================================
addpath("Utilities/");
addpath("Calibration/");
addpath("ex2/");
addpath("Distributions/");

% Run the full calibration (verbose report + implied-distribution plot).
opts = struct();
opts.valueDate = datetime(2020, 06, 02);
opts.verbose   = true;
opts.plot      = true;

[params, market] = calibrate_surface(opts);

% Calibrated parameters are now available in `params`:
%   params.AB.k, params.AB.eta, params.AB.sigma_t, ...
%   params.MA.alpha, params.MA.beta, ...
%   params.GL.alpha, params.GL.beta, ...
% and the market/support data (forward, discount_factor, sigma_ATM, yf, ...)
% in `market`.
