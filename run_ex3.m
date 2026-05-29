%% MAIN SCRIPT - EX 3: FORWARD START PRICING  (MA + AB + GL)
% Copy of run_exercise_3 (Minimal Additive) extended with the Additive
% Bachelier (AB) and Generalized Logistic (GL) models. The market curve, ATM
% volatility and the forward-start window [t1 = yf(2), t2 = yf(4)] are shared;
% each model is then priced with its OWN native functions:
%   MA -> Functions_ex3   (FA_simulation, pricing_fwd_start_*)
%   AB -> ex2 + ex3_utilities (calibrateAB, ccdf_AB_FFT, price_AB_MC)
%   GL -> ex3_GL          (lewis_fft_cdf_old, simulate_increments, cf_increment_GL)
%
% For each model we compare the forward-start price (analytic / FFT reference)
% against Monte Carlo. MA is checked across many strikes; AB and GL at the
% ATM forward start (K2 = 1), where the payoff reduces to max(W, 0).

clear; clc; close all;

% =========================================================================
% PATHS INITIALIZATION
% =========================================================================
addpath("Utilities/");
addpath("Functions/");          % I0_MA, cf_GL, pdf_GL, complex_gamma
addpath("Functions_ex3/");      % MA engine
addpath("ex2/");                % calibrateAB, call_AB_FFT, I0, charateristic_function_AB
addpath("ex3_utilities/");      % ccdf_AB_FFT, sample_from_cdf, price_AB_MC, plot_mc_check
addpath("ex3_GL/");             % lewis_fft_cdf_old, simulate_increments, cf_increment_GL

fprintf('=========================================================================\n');
fprintf('        LINEAR ADDITIVE MODELS (MA / AB / GL) - FORWARD-START ENGINE      \n');
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

% Print Bootstrap Results
fprintf("\n--- Bootstrap Results (Value Date: %s) ---\n", string(valueDate, "yyyy-MM-dd"));
fprintf("%-12s  %10s  %10s  %8s\n", "Expiry", "D(t,T)", "F(t,T)", "R^2");
for k = 1:nT
    fprintf("%-12s  %10.6f  %10.4f  %8.4f\n", ...
        string(expiries(k), "yyyy-MM-dd"), ...
        discount_factor(k), forward(k), R2(k));
end
fprintf('-------------------------------------------------------------------------\n\n');

%% =========================================================================
% STEP 2: ATM VOLATILITY CALIBRATION
% =========================================================================
fprintf('STEP 2: Calibrating ATM Volatility...\n');
c_ATM = zeros(length(forward), 1);

for i = 1:length(forward)
    current_calls = calls(i, :);
    current_puts  = puts(i, :);
    c_ATM(i)      = callATM(current_calls, current_puts, strikes, forward(i), discount_factor(i));
end

yf = yearfrac(valueDate, expiries, 3);
sigma_ATM = sigmaATM(c_ATM, discount_factor, yf, expiries);
fprintf('  -> ATM Volatility successfully calibrated for %d maturities.\n\n', nT);

% Shared forward-start window: t1 = yf(2), t2 = yf(4)
iT1 = 2;  iT2 = 4;
T1  = yf(iT1);  T2 = yf(iT2);
fprintf('  Forward-start window: t1 = %.4f y (idx %d), t2 = %.4f y (idx %d)\n\n', T1, iT1, T2, iT2);

%% =========================================================================
% STEP 3: SANITY CHECK - FFT VS ANALYTICAL GAUSSIAN CDF
% =========================================================================
fprintf('STEP 3: Validating FFT numerical inversion against Gaussian CDF...\n');
M_gauss     = 14;
dz_gauss    = 2.5e-3;
shift_gauss = 0.5;

[cdf_clean_gauss, z_grid_gauss] = lewis_fft_algorithm_trial(M_gauss, dz_gauss, shift_gauss, 1);
cdf_analytical = normcdf(z_grid_gauss);

figure('Name', 'FFT vs Analytical Gaussian CDF');
plot(z_grid_gauss, cdf_clean_gauss, 'b-', 'LineWidth', 2, 'DisplayName', 'FFT Inversion');
hold on;
plot(z_grid_gauss, cdf_analytical, 'r--', 'LineWidth', 2, 'DisplayName', 'Analytical CDF');
title('FFT Inversion vs Analytical Gaussian CDF');
xlabel('Spatial Grid (z)'); ylabel('Cumulative Probability');
legend('Location', 'northwest'); xlim([-5, 5]); grid on; hold off;
fprintf('  -> Validation complete.\n\n');

%% #########################################################################
%  MODEL 1 - MINIMAL ADDITIVE (MA)
%  #########################################################################

%% =========================================================================
% STEP 4 (MA): MODEL PARAMETERS SETUP
% =========================================================================
fprintf('STEP 4 (MA): Setting up Minimal Additive parameters...\n');
rng(2); % For reproducibility

% 1. Model Parameters definition (calibrated in Exercise 2)
alpha_MA = 4.90;
beta_MA  = 5.28;
gamma_MA = (1/alpha_MA) - (1/beta_MA);

% 2. Relevant dates / vols on the forward-start window
yf_fwd        = [yf(iT1); yf(iT2)];
sigma_ATM_fwd = [sigma_ATM(iT1); sigma_ATM(iT2)];
df_fwd        = discount_factor(iT2);
forward_fwd   = forward(iT2);

% 3. Integrated Volatility computation
% NOTE: I0_MA has signature I0_MA(gamma, C, alpha, beta) with C = (1/a+1/b)^-1
% (see price_MA.m). The original run_exercise_3 called it with 2 args -> fixed.
C_MA   = 1 / ((1/alpha_MA) + (1/beta_MA));
I0     = I0_MA(alpha_MA, beta_MA);
sigmat = (sigma_ATM_fwd / I0) .* yf_fwd;

fprintf('  -> Parameters extracted for interval [t1 = %.2f, t2 = %.2f]\n\n', yf_fwd(1), yf_fwd(2));

%% =========================================================================
% STEP 5 (MA): FORWARD START OPTION PRICING (MULTI-STRIKE)
% =========================================================================
fprintf('STEP 5 (MA): Pricing Forward Start Option (Analytic vs Monte Carlo)...\n');

N_sim_MA = 1e7;
M_MA     = 16;
dz_MA    = 5e-3;

% Proportional strikes: ATM (1.0) + 20 random in [0.8, 1.2]
num_random_strikes = 20;
K2_vec = [1.0, 0.8 + 0.4 * rand(1, num_random_strikes)];
K2_vec = sort(K2_vec);

% Analytical Pricing (vectorized)
price_analytic = pricing_fwd_start_analytic(alpha_MA, beta_MA, sigmat, df_fwd, K2_vec, forward_fwd);

% Monte Carlo Pricing (vectorized)
fprintf('  -> Running Monte Carlo Simulation (%d paths)...\n', N_sim_MA);
[price_MC, CI_MC] = pricing_fwd_start_MC(forward_fwd, K2_vec, df_fwd, N_sim_MA, M_MA, dz_MA, sigmat, alpha_MA, beta_MA);

% Comparison table
diff_bps = (price_MC - price_analytic) * 10000;
T_MA = table(K2_vec', price_analytic', price_MC', CI_MC(1,:)', CI_MC(2,:)', diff_bps', ...
    'VariableNames', {'Strike_Multiplier', 'Analytic_Price', 'MC_Price', 'CI_Lower', 'CI_Upper', 'Diff_bps'});
fprintf('\n--- MA FORWARD START: ANALYTIC VS MONTE CARLO ---\n');
disp(T_MA);

% Plot
idx_atm = find(K2_vec == 1);
figure('Name', 'MA Forward Start Pricing');
hold on; grid on;
plot(K2_vec, price_analytic, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytic Price');
plot(K2_vec, price_MC, 'ro', 'MarkerSize', 6, 'DisplayName', 'Monte Carlo Price');
if ~isempty(idx_atm)
    plot(K2_vec(idx_atm), price_analytic(idx_atm), 'kp', ...
        'MarkerSize', 15, 'MarkerFaceColor', 'y', 'DisplayName', 'ATM Forward (K_2 = 1)');
end
title('MA Forward Start Pricing: Analytic vs Monte Carlo');
xlabel('Proportional Strike (K_2)'); ylabel('Discounted Option Price');
legend('Location', 'best'); xlim([min(K2_vec), max(K2_vec)]); hold off;

%% =========================================================================
% STEP 6 (MA): ROBUSTNESS CHECK - EXACT CDF NUMERICAL INTEGRATION (ATM)
% =========================================================================
fprintf('\nSTEP 6 (MA): Verifying ATM Analytical Price via Numerical Integration...\n');
if ~isempty(idx_atm)
    price_analytic_atm = price_analytic(idx_atm);

    ps_plus  = beta_MA / sigmat(1);
    ps_minus = alpha_MA / sigmat(1);
    pt_plus  = beta_MA / sigmat(2);
    pt_minus = alpha_MA / sigmat(2);
    drift_t1_t2 = gamma_MA * (sigmat(2) - sigmat(1));

    [cdf_exact, x_grid] = exact_ma_increment_cdf(pt_plus, pt_minus, ps_plus, ps_minus, drift_t1_t2);
    cdf_fun = @(x) interp1(x_grid, cdf_exact, x, 'pchip', 'extrap');
    expectation_analytical = integral(@(x) 1 - cdf_fun(x), 0, x_grid(end));

    price_integrated = df_fwd * expectation_analytical;
    diff_bps_integration = abs(price_integrated - price_analytic_atm) * 10000;

    fprintf('  -> Formula-based Price (K2=1) : %.6f\n', price_analytic_atm);
    fprintf('  -> Integral-based Price       : %.6f\n', price_integrated);
    fprintf('  -> Integration Error          : %.2f bps\n\n', diff_bps_integration);
end

%% #########################################################################
%  MODEL 2 - ADDITIVE BACHELIER (AB)
%  #########################################################################
fprintf('STEP 7 (AB): Calibrating (k, eta) and pricing forward start...\n');

kAB=0.93;
eta_AB=-0.06;
I_0=I0_AB(0,kAB,eta_AB);
sigma_t_AB = sigma_ATM / I_0; % Rescaling the forward vols by the MA normalization I0

% --- conditional CDF of the increment t1 -> t2 (Lemma 2 forward rescaling) ---
N_sim_AB   = 1e6;
sigma_T1   = sigma_t_AB(iT1);   sigma_T2 = sigma_t_AB(iT2);
fwd_factor = discount_factor(iT1) / discount_factor(iT2);

% grid scaled to the increment std: Var(W) = (1+k*eta^2)*(sT2^2*T2 - fwd^2*sT1^2*T1)
std_W      = sqrt((1 + kAB*eta_AB^2) * (sigma_T2^2*T2 - fwd_factor^2*sigma_T1^2*T1));
x_grid_AB  = linspace(-10*std_W, 10*std_W, 2000)';
cdf_AB     = ccdf_AB_FFT(eta_AB, kAB, T1, T2, sigma_T1, sigma_T2, x_grid_AB, 0, fwd_factor);

% --- MC simulation + visual check ---
Z_AB = sample_from_cdf(x_grid_AB, cdf_AB, N_sim_AB);
plot_mc_check(Z_AB, x_grid_AB, cdf_AB, T1, T2);

% --- forward-start K2 = 1: MC vs FFT (survival-integral) reference ---
[price_AB_mc, IC_AB] = price_AB_MC(T1, T2, kAB, eta_AB, sigma_T1, sigma_T2, N_sim_AB, ...
                                   forward(iT2), discount_factor(iT1), discount_factor(iT2), x_grid_AB, 1);

F0_AB    = interp1(x_grid_AB, cdf_AB, 0, 'spline');   % F_W(0)
mask_pos = x_grid_AB > 0;
price_AB_fft = discount_factor(iT2) * trapz([0; x_grid_AB(mask_pos)], 1 - [F0_AB; cdf_AB(mask_pos)]);

fprintf('\n--- AB FORWARD START (K2 = 1) ---\n');
fprintf('  FFT reference (analytic) : %.6f\n', price_AB_fft);
fprintf('  Monte Carlo              : %.6f   95%% CI [%.6f, %.6f]\n', price_AB_mc, IC_AB(1), IC_AB(2));
fprintf('  Difference               : %.2f bps\n\n', abs(price_AB_mc - price_AB_fft)*1e4);

%% #########################################################################
%  MODEL 3 - GENERALIZED LOGISTIC (GL)
%  #########################################################################
fprintf('STEP 8 (GL): Pricing forward start...\n');

% Model parameters (calibrated in Exercise 2)
alpha_GL = 0.44;
beta_GL  = 0.40;

% GL normalization I0 = sqrt(2*pi)*E[zeta_+]  (Baviera-Massaria Eq.14)
integrand_mean = @(x) pdf_GL(alpha_GL, beta_GL, x) .* x;
I0_GL = sqrt(2*pi) * quadgk(integrand_mean, 0, inf);

% Scales on the forward-start window
sigma_s_GL = sigma_ATM(iT1) / I0_GL;
sigma_t_GL = sigma_ATM(iT2) / I0_GL;

% Conditional CDF of the increment t1 -> t2 via Lewis-FFT (double-shift splicing)
[CDF_GL, zk_GL] = lewis_fft_cdf_old(@cf_increment_GL, alpha_GL, beta_GL, ...
                                    sigma_s_GL, T1, sigma_t_GL, T2);

% Simulate the t1 -> t2 increment (inverse-CDF)
N_sim_GL = 1e6;  seed_GL = 42;
[X_st_GL, x_grid_GL, CDF_clean_GL] = simulate_increments(zk_GL, CDF_GL, N_sim_GL, seed_GL);

% Forward-start K2 = 1: payoff = max(S_t - F(t1,t2), 0) = max(W, 0)
price_GL_mc = discount_factor(iT2) * mean(max(X_st_GL, 0));

% Analytic (FFT) reference: E[max(W,0)] = integral_0^inf (1 - F_W) dx
cdf_fun_GL  = @(x) interp1(x_grid_GL, CDF_clean_GL, x, 'spline', 'extrap');
price_GL_an = discount_factor(iT2) * integral(@(x) 1 - cdf_fun_GL(x), 0, x_grid_GL(end));

fprintf('\n--- GL FORWARD START (K2 = 1) ---\n');
fprintf('  FFT reference (analytic) : %.6f\n', price_GL_an);
fprintf('  Monte Carlo              : %.6f\n', price_GL_mc);
fprintf('  Difference               : %.2f bps\n\n', abs(price_GL_mc - price_GL_an)*1e4);

%% =========================================================================
% FINAL SUMMARY (K2 = 1)
% =========================================================================
fprintf('=========================================================================\n');
fprintf('              FORWARD-START SUMMARY AT K2 = 1 (ATM)                       \n');
fprintf('=========================================================================\n');
if ~isempty(idx_atm)
    fprintf('  MA : analytic %.6f | MC %.6f\n', price_analytic(idx_atm), price_MC(idx_atm));
end
fprintf('  AB : analytic %.6f | MC %.6f\n', price_AB_fft, price_AB_mc);
fprintf('  GL : analytic %.6f | MC %.6f\n', price_GL_an, price_GL_mc);
fprintf('=========================================================================\n');
fprintf('                 SIMULATION COMPLETED SUCCESSFULLY.                       \n');
fprintf('=========================================================================\n');
