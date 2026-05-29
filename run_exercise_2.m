%% MAIN SCRIPT - VOLATILITY SURFACE CALIBRATION (MA & GL MODELS)
% This script loads market option data, bootstraps the forward curve and 
% discount factors, calibrates the ATM volatility, and fits both the 
% Minimal Additive (MA) and Generalized Laplace (GL) models to the market surface.

clear; clc; close all;

% =========================================================================
% PATHS INITIALIZATION
% =========================================================================
addpath("Utilities/");
addpath("Calibration/");

fprintf('=========================================================================\n');
fprintf('             VOLATILITY SURFACE CALIBRATION ENGINE                       \n');
fprintf('=========================================================================\n\n');

%% =========================================================================
% STEP 1: DATA LOADING AND BOOTSTRAPPING
% =========================================================================
fprintf('STEP 1: Loading data and bootstrapping curve...\n');
callpath   = "Data/datacalls";
putpath    = "Data/dataputs";
expiryFile = "Data/Expiries_Futures.txt";
valueDate  = datetime(2020, 06, 02);

% Read market options data
[strikes, calls, puts, expiries] = readData(callpath, putpath, valueDate, expiryFile);

% Bootstrap synthetic discount factors and forwards
nT = numel(expiries);
discount_factor = zeros(nT, 1);
forward         = zeros(nT, 1);
R2              = zeros(nT, 1);

for k = 1:nT
    [discount_factor(k), forward(k), R2(k)] = bootstrap(puts(k,:), calls(k,:), strikes);
end

% Report Bootstrap Results
fprintf("\n--- Bootstrap Results (Value Date: %s) ---\n", string(valueDate, "yyyy-MM-dd"));
fprintf("%-12s  %10s  %10s  %8s\n", "Expiry", "D(t,T)", "F(t,T)", "R^2");
for k = 1:nT
    fprintf("%-12s  %10.6f  %10.4f  %8.4f\n", ...
        string(expiries(k), "yyyy-MM-dd"), ...
        discount_factor(k), forward(k), R2(k));
end
fprintf('-------------------------------------------------------------------------\n\n');

%% =========================================================================
% STEP 2: ATM VOLATILITY CALIBRATION & MONEYNESS GENERATION
% =========================================================================
fprintf('STEP 2: Calibrating ATM Volatility and filtering surface...\n');
c_ATM = zeros(length(forward), 1); 

for i = 1:length(forward)
    % Extract data for current maturity
    current_calls = calls(i, :);
    current_puts  = puts(i, :);
    
    c_ATM(i) = callATM(current_calls, current_puts, strikes, forward(i), discount_factor(i));
end

% Since expiries are already business-day adjusted, we calculate the standard year fraction
yf = yearfrac(valueDate, expiries, 3);
sigma_ATM = sigmaATM(c_ATM, discount_factor, yf, expiries);

% Verify the additive-model precondition (sigma_ATM * sqrt(t) increasing).
% No (alpha, beta) calibration can fix this if it fails -- it is a data
% condition. The check only warns; calibration still proceeds.
check_term_structure(sigma_ATM, yf, expiries);

% Generate the modified moneyness matrix and filter out-of-bounds market prices.
% Bounds are now in *normalized* moneyness chi = (K - F)/(sigma_ATM*sqrt(t)),
% so the band covers a comparable number of standard deviations at every
% maturity (cf. paper 3 Sec. 3, separability of the implied vol).
chi_min = -5;
chi_max = 5;
[moneyness_modified, c_mkt_calibration] = moneyness_generator( ...
    forward, strikes, calls, puts, sigma_ATM, yf, discount_factor, ...
    chi_min, chi_max);
fprintf('  -> Surface filtered on chi in [%g, %g] and prepared for optimization.\n\n', ...
    chi_min, chi_max);

%% =========================================================================
% STEP 3: ADDITIVE BACHELIER (AB) MODEL CALIBRATION
% =========================================================================
fprintf('STEP 3: Calibrating Additive Bachelier (AB) Model via fmincon...\n');

% AB is now calibrated exactly like MA and GL: directly on the (M x N)
% matrices returned by moneyness_generator (moneyness_modified,
% c_mkt_calibration), with NaNs marking strikes outside the band / without a
% market quote. No per-maturity for-loop and no pooled vectors -- price_AB
% prices the whole surface in a single Lewis-FFT (the normalized AB price
% G(chi; eta, k) is maturity-independent under separability, paper 3 Eq. 19)
% and the dollar price is recovered row-wise via B * sigma_ATM * sqrt(t)
% (paper 3 Eq. 20, the same loss used in run_project2A).
%
% Parameters x = [k, eta]. Bounds from paper 3 Fig. 2-3 (Covid sample):
% eta ~ [-0.3, 0.3], k ~ [0.4, 1.2]; widened with margin and with k kept away
% from 0, where the Lewis contour strip collapses.
options_AB = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'interior-point', ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance',       1e-10, ...
    'MaxFunctionEvaluations', 5000);

x0_AB = [1.0, 0.2];
lb_AB = [1e-3, -1.5];
ub_AB = [5.0,   1.5];

obj_fun_AB = @(x) objective_function_AB(x, discount_factor, yf, sigma_ATM, ...
                       moneyness_modified, c_mkt_calibration);

[x_opt_AB, fval_AB, exitflag_AB] = ...
    fmincon(obj_fun_AB, x0_AB, [], [], [], [], lb_AB, ub_AB, [], options_AB);

k_AB   = x_opt_AB(1);
eta_AB = x_opt_AB(2);

fprintf('\n  -> Optimal parameters found: k = %.6f, eta = %.6f.\n', k_AB, eta_AB);
fprintf('  -> AB Calibration completed (exitflag = %d, SSE = %.6g).\n\n', ...
        exitflag_AB, fval_AB);

%% =========================================================================
% STEP 3: MINIMAL ADDITIVE (MA) MODEL CALIBRATION
% =========================================================================
fprintf('STEP 3: Calibrating Minimal Additive (MA) Model via fminbnd...\n');

%% MA: ft=sigma_ATM * sqrt(t) * zeta_t dove zeta_t asymmetric laplace distribution

% --- MA scale invariance ---------------------------------------------------
% The MA call price is exactly invariant under (alpha, beta) -> (s*alpha,
% s*beta) for any s > 0: gamma -> gamma/s, C -> s*C, I0 -> I0/s, so
% maturity_multiplier scales by s while `core` scales by 1/s and the two
% cancel. Hence the SSE landscape has a perfectly flat 1-D ridge along
% every line through the origin, and only the *ratio* beta/alpha is
% identifiable from market prices. To remove the redundancy we fix
% alpha = 1 and calibrate beta in a 1-D search; the overall scale is then
% absorbed by the cascade definition sigma_MA = sigma_ATM / I0.
ALPHA_FIX = 1.0;

lb_beta = 0.05;
ub_beta = 40;

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'interior-point', ...
    'FunctionTolerance',   1e-12, ...
    'StepTolerance',       1e-12, ...
    'OptimalityTolerance', 1e-10, ...
    'MaxFunctionEvaluations', 2000);

obj_fun_MA_1d = @(b) objective_function_MA([ALPHA_FIX, b], ...
                       discount_factor, yf, sigma_ATM, ...
                       moneyness_modified, c_mkt_calibration);

% MA is a smooth, 1-D, bounded problem (alpha fixed as gauge, only beta is
% free, since the scale invariance makes only beta/alpha identifiable).
% fminbnd is purpose-built for it: golden-section + parabolic interpolation,
% no initial guess, stays inside [lb_beta, ub_beta] by construction, and
% avoids fmincon's finite-difference gradients near the penalty cliff in the
% objective. The term-structure nonlcon is omitted on purpose: with
% (alpha, beta) constant it does not depend on beta, so it cannot bind (the
% data condition is already reported by check_term_structure above).
opts_MA_1d = optimset('TolX', 1e-10, 'Display', 'iter');

[beta_opt, fval_MA, exitflag_MA] = ...
    fminbnd(obj_fun_MA_1d, lb_beta, ub_beta, opts_MA_1d);

alpha_MA = ALPHA_FIX;
beta_MA  = beta_opt;

fprintf('\n  -> MA Calibration completed (exitflag = %d, SSE = %.6g).\n', ...
        exitflag_MA, fval_MA);
fprintf('  -> Fixed scale: alpha = %.6f (gauge); calibrated beta = %.6f.\n', ...
        alpha_MA, beta_MA);
fprintf('  -> Asymmetry ratio beta/alpha = %.6f.\n\n', beta_MA/alpha_MA);

% --- MA calibration check via implied-vol SKEW -----------------------------
% The asymmetric-Laplace skewness is scale-invariant and depends only on
% beta/alpha (the single identifiable MA parameter), and under separability
% the MA smile in chi is maturity-flat. So matching the market skew is the
% most direct sanity check on the MA fit. See check_skew_MA.m.
skew_report_MA = check_skew_MA(alpha_MA, beta_MA, discount_factor, yf, ...
    sigma_ATM, moneyness_modified, c_mkt_calibration, expiries);

% Re-define nonlcon for GL (2-D) -- it does not have the MA scale invariance,
% so it is genuinely calibrated on both alpha and beta.
nonlcon_TS = @(x) term_structure_nonlcon(x, sigma_ATM, yf);

%% =========================================================================
% STEP 4: GENERALIZED LAPLACE (GL) MODEL CALIBRATION
% =========================================================================
fprintf('STEP 4: Calibrating Generalized Laplace (GL) Model via fmincon...\n');

% FFT Numerical settings
M  = 15;
dz = 2.5e-3;

% Optimization settings for GL
x0_GL = [0.9, 0.5];              
lb_GL = [0.05, 0.05];             
ub_GL = [50,  50];                

obj_fun_GL = @(x) objective_function_GL(x, discount_factor, yf, sigma_ATM, moneyness_modified, c_mkt_calibration, M, dz);
[x_opt_GL, fval_GL, exitflag_GL] = fmincon(obj_fun_GL, x0_GL, [], [], [], [], lb_GL, ub_GL, nonlcon_TS, options);

alpha_GL = x_opt_GL(1);
beta_GL  = x_opt_GL(2);
fprintf('\n  -> Optimal parameters found: alpha = %.6f, beta = %.6f.\n', alpha_GL, beta_GL);
fprintf('\n  -> GL Calibration completed (exitflag = %d, SSE = %.6g).\n\n', exitflag_GL, fval_GL);

%% =========================================================================
% STEP 5: POST-CALIBRATION DIAGNOSTICS & SUMMARY
% =========================================================================
fprintf('STEP 5: Post-Calibration Diagnostics and Reporting...\n');

% Compute model prices using calibrated parameters
c_mod_AB = price_AB([k_AB, eta_AB], discount_factor, yf, sigma_ATM, moneyness_modified);
c_mod_MA = price_MA([alpha_MA, beta_MA], discount_factor, yf, sigma_ATM, moneyness_modified);
c_mod_GL = price_GL(alpha_GL, beta_GL, M, dz, discount_factor, sigma_ATM, yf, moneyness_modified);

% Clamp tiny negative artefacts from FFT/interp extrapolation to 0.
% (NaN < 0 is false in MATLAB, so NaN entries are preserved.)
c_mod_AB(c_mod_AB < 0) = 0;
c_mod_MA(c_mod_MA < 0) = 0;
c_mod_GL(c_mod_GL < 0) = 0;

% Calculate residuals
res_AB = c_mkt_calibration - c_mod_AB;
res_MA = c_mkt_calibration - c_mod_MA;
res_GL = c_mkt_calibration - c_mod_GL;

% Print Parameter Summary Table
fprintf('\n=========================================================================\n');
fprintf('                     CALIBRATED PARAMETERS SUMMARY                       \n');
fprintf('=========================================================================\n');
fprintf('%-10s | %-12s | %-12s | %-15s\n', 'Model', 'Param 1', 'Param 2', 'Final SSE');
fprintf('-------------------------------------------------------------------------\n');
fprintf('%-10s | %-12.6f | %-12.6f | %-15.6g\n', 'AB (k,eta)', k_AB, eta_AB, fval_AB);
fprintf('%-10s | %-12.6f | %-12.6f | %-15.6g\n', 'MA (a,b)', alpha_MA, beta_MA, fval_MA);
fprintf('%-10s | %-12.6f | %-12.6f | %-15.6g\n', 'GL (a,b)', alpha_GL, beta_GL, fval_GL);

% Print RMSE per Maturity Table
fprintf('\n=========================================================================\n');
fprintf('                     ROOT MEAN SQUARE ERROR (RMSE)                       \n');
fprintf('=========================================================================\n');
fprintf('%-15s | %-12s | %-12s | %-12s\n', 'Expiry Date', 'RMSE AB', 'RMSE MA', 'RMSE GL');
fprintf('-------------------------------------------------------------------------\n');
for i = 1:nT
    r_AB_i = res_AB(i,:); r_AB_i = r_AB_i(~isnan(r_AB_i));
    r_MA_i = res_MA(i,:); r_MA_i = r_MA_i(~isnan(r_MA_i));
    r_GL_i = res_GL(i,:); r_GL_i = r_GL_i(~isnan(r_GL_i));
    fprintf('%-15s | %-12.6f | %-12.6f | %-12.6f\n', string(expiries(i),'yyyy-MM-dd'), ...
        sqrt(mean(r_AB_i.^2)), sqrt(mean(r_MA_i.^2)), sqrt(mean(r_GL_i.^2)));
end

% Print Sample Pricing Comparison Table (Extracting ~10 valid non-NaN prices)
fprintf('\n=============================================================================================================\n');
fprintf('                             SAMPLE PRICING COMPARISON (MARKET VS MODELS)                                    \n');
fprintf('=============================================================================================================\n');
fprintf('%-10s | %-10s | %-12s | %-12s | %-12s | %-12s | %-13s | %-13s | %-13s\n', 'Maturity', 'Moneyness', 'Market Px', 'AB Price', 'MA Price', 'GL Price', 'Diff AB (bps)', 'Diff MA (bps)', 'Diff GL (bps)');
fprintf('-------------------------------------------------------------------------------------------------------------\n');

% Find valid indices to display a subset of results
[mat_idx, mon_idx] = find(~isnan(c_mkt_calibration));
num_samples = min(10, length(mat_idx)); 

% Select roughly evenly spaced samples from the valid indices pool
sample_indices = round(linspace(1, length(mat_idx), num_samples));

for k = 1:num_samples
    idx = sample_indices(k);
    r = mat_idx(idx);
    c = mon_idx(idx);
    
    % Calculate differences in basis points (Model - Market)
    diff_AB_bps = (c_mod_AB(r, c) - c_mkt_calibration(r, c)) * 10000;
    diff_MA_bps = (c_mod_MA(r, c) - c_mkt_calibration(r, c)) * 10000;
    diff_GL_bps = (c_mod_GL(r, c) - c_mkt_calibration(r, c)) * 10000;

    fprintf('%-10s | %-10.4f | %-12.6f | %-12.6f | %-12.6f | %-12.6f | %-13.2f | %-13.2f | %-13.2f\n', ...
        string(expiries(r), 'yyyy-MM'), ...
        moneyness_modified(r, c), ...
        c_mkt_calibration(r, c), ...
        c_mod_AB(r, c), ...
        c_mod_MA(r, c), ...
        c_mod_GL(r, c), ...
        diff_AB_bps, ...
        diff_MA_bps, ...
        diff_GL_bps);
end
fprintf('=============================================================================================================\n');
%% =========================================================================
% STEP 6: PLOTTING IMPLED DISTRIBUTIONS
% =========================================================================
% Use the AB parameters calibrated above (k_AB is the IG subordinator
% parameter kappa in plot_distributions' NIG mapping).
plot_distributions(eta_AB, k_AB, alpha_MA, beta_MA, alpha_GL, beta_GL);