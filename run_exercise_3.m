%% MAIN SCRIPT - MINIMAL ADDITIVE (MA) PRICING & SIMULATION
% This script loads market data, calibrates the ATM volatility, validates
% the FFT numerical inversion, and prices a Forward Start Option using 
% both Monte Carlo simulation and exact analytical formulas under the 
% Minimal Additive (MA) framework across multiple strikes.

clear; clc; close all;

% =========================================================================
% PATHS INITIALIZATION
% =========================================================================
addpath("Utilities/");
addpath("Functions_ex2/");
addpath("Functions_ex3/");

fprintf('=========================================================================\n');
fprintf('             MINIMAL ADDITIVE (MA) MODEL - PRICING ENGINE                \n');
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

% Optional: Generate moneyness and market prices for further calibration
[moneyness_mod, c_mkt_calib] = moneyness_generator(forward, strikes, calls, puts, sigma_ATM, yf, discount_factor, -30, 30);
fprintf('  -> ATM Volatility successfully calibrated for %d maturities.\n\n', nT);

%% =========================================================================
% STEP 3: SANITY CHECK - FFT VS ANALYTICAL GAUSSIAN CDF
% =========================================================================
fprintf('STEP 3: Validating FFT numerical inversion against Gaussian CDF...\n');
M_gauss     = 14;
dz_gauss    = 2.5e-3;
shift_gauss = 0.5;

% 1. Compute Numerical CDF via FFT
[cdf_clean_gauss, z_grid_gauss] = lewis_fft_algorithm_trial(M_gauss, dz_gauss, shift_gauss, 1);

% 2. Compute Exact Analytical CDF
cdf_analytical = normcdf(z_grid_gauss);

% 3. Plot Comparison
figure('Name', 'FFT vs Analytical Gaussian CDF');
plot(z_grid_gauss, cdf_clean_gauss, 'b-', 'LineWidth', 2, 'DisplayName', 'FFT Inversion');
hold on;
plot(z_grid_gauss, cdf_analytical, 'r--', 'LineWidth', 2, 'DisplayName', 'Analytical CDF');
title('FFT Inversion vs Analytical Gaussian CDF');
xlabel('Spatial Grid (z)'); ylabel('Cumulative Probability');
legend('Location', 'northwest'); xlim([-5, 5]); grid on; hold off;
fprintf('  -> Validation complete. Check the generated plot.\n\n');

%% =========================================================================
% STEP 4: MODEL PARAMETERS SETUP
% =========================================================================
fprintf('STEP 4: Setting up Minimal Additive (MA) Model Parameters...\n');
rng(2); % For reproducibility

% 1. Model Parameters definition
alpha_MA = 4.90;
beta_MA  = 5.28;
gamma_MA = (1/alpha_MA) - (1/beta_MA);

% 2. Extraction of relevant dates (t0=0, t1=yf(2), t2=yf(4))
yf_fwd        = [yf(2); yf(4)];
sigma_ATM_fwd = [sigma_ATM(2); sigma_ATM(4)];
df_fwd        = discount_factor(4);
forward_fwd   = forward(4);

% 3. Integrated Volatility computation
I0 = I0_MA(alpha_MA, beta_MA);
sigmat = (sigma_ATM_fwd / I0) .* yf_fwd;

fprintf('  -> Parameters extracted for interval [t1 = %.2f, t2 = %.2f]\n\n', yf_fwd(1), yf_fwd(2));

%% =========================================================================
% STEP 5: FORWARD START OPTION PRICING (MULTI-STRIKE COMPARISON)
% =========================================================================
fprintf('STEP 5: Pricing Forward Start Option (Analytic vs Monte Carlo)...\n');

% 1. Grid & Simulation settings
N_sim = 1e7;
M     = 16;
dz    = 5e-3;

% 2. Define Proportional Strikes (K2)
% We include the ATM-Forward case (1.0) and add 20 random strikes in [0.8, 1.2]
num_random_strikes = 20;
K2_vec = [1.0, 0.8 + 0.4 * rand(1, num_random_strikes)];
K2_vec = sort(K2_vec); % Sort for clean tabular display and plotting

% 3. Analytical Pricing (Vectorized)
price_analytic = pricing_fwd_start_analytic(alpha_MA, beta_MA, sigmat, df_fwd, K2_vec, forward_fwd);

% 4. Monte Carlo Pricing (Vectorized)
fprintf('  -> Running Monte Carlo Simulation (%d paths)...\n', N_sim);
[price_MC, CI_MC] = pricing_fwd_start_MC(forward_fwd, K2_vec, df_fwd, N_sim, M, dz, sigmat, alpha_MA, beta_MA);

% 5. Print Comparison Table
diff_bps = (price_MC - price_analytic) * 10000;
T_comparison = table(K2_vec', price_analytic', price_MC', CI_MC(1,:)', CI_MC(2,:)', diff_bps', ...
    'VariableNames', {'Strike_Multiplier', 'Analytic_Price', 'MC_Price', 'CI_Lower', 'CI_Upper', 'Diff_bps'});

fprintf('\n=========================================================================\n');
fprintf('       FORWARD START COMPARISON: ANALYTIC VS MONTE CARLO                 \n');
fprintf('=========================================================================\n');
disp(T_comparison);

% 6. Plotting the Results
figure('Name', 'Forward Start Pricing');
hold on; grid on;

% Plot Analytical Curve
plot(K2_vec, price_analytic, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytic Price');

% Plot Monte Carlo Points
plot(K2_vec, price_MC, 'ro', 'MarkerSize', 6, 'DisplayName', 'Monte Carlo Price');

% Highlight the ATM Forward (K2 = 1.0) with a star
idx_atm = find(K2_vec == 1);
if ~isempty(idx_atm)
    plot(K2_vec(idx_atm), price_analytic(idx_atm), 'kp', ...
        'MarkerSize', 15, 'MarkerFaceColor', 'y', 'DisplayName', 'ATM Forward (K_2 = 1)');
end

title('Forward Start Option Pricing: Analytic vs Monte Carlo', 'FontSize', 12);
xlabel('Proportional Strike (K_2)', 'FontSize', 11);
ylabel('Discounted Option Price', 'FontSize', 11);
legend('Location', 'best');
xlim([min(K2_vec), max(K2_vec)]);
hold off;

%% =========================================================================
% STEP 6: ROBUSTNESS CHECK - EXACT CDF NUMERICAL INTEGRATION (ATM CASE)
% =========================================================================
fprintf('\nSTEP 6: Verifying ATM Analytical Price via Numerical Integration...\n');

% Check if the ATM price is available in our computed vectors
if ~isempty(idx_atm)
    price_analytic_atm = price_analytic(idx_atm);
    
    % Tail decays and drift for the increment [t1, t2]
    ps_plus  = beta_MA / sigmat(1);
    ps_minus = alpha_MA / sigmat(1);
    pt_plus  = beta_MA / sigmat(2);
    pt_minus = alpha_MA / sigmat(2);
    drift_t1_t2 = gamma_MA * (sigmat(2) - sigmat(1));

    % 1. Get the exact analytical CDF over a spatial grid
    [cdf_exact, x_grid] = exact_ma_increment_cdf(pt_plus, pt_minus, ps_plus, ps_minus, drift_t1_t2);

    % 2. Create a continuous function handle using spline extrapolation
    cdf_fun = @(x) interp1(x_grid, cdf_exact, x, 'pchip', 'extrap'); 

    % 3. Integrate the survival function: E[X^+] = \int_0^\infty (1 - F(x)) dx
    integrand = @(x) 1 - cdf_fun(x);
    expectation_analytical = integral(integrand, 0, x_grid(end));

    % 4. Discount back to present
    price_integrated = df_fwd * expectation_analytical;
    diff_bps_integration = abs(price_integrated - price_analytic_atm) * 10000;

    % Print Integration Results
    fprintf('  -> Formula-based Price (K2=1) : %.6f\n', price_analytic_atm);
    fprintf('  -> Integral-based Price       : %.6f\n', price_integrated);
    fprintf('  -> Integration Error          : %.2f bps\n', diff_bps_integration);
else
    fprintf('  -> SKIPPED: ATM Strike (K2=1) not found in the evaluation vector.\n');
end

fprintf('=========================================================================\n');
fprintf('                 SIMULATION COMPLETED SUCCESSFULLY.                      \n');
fprintf('=========================================================================\n');