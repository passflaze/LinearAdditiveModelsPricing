function LA_results = run_ex3(params, market)
%% EX 3: FORWARD START PRICING  (MA + AB + GL)
% run_exercise_3 for the 3 models, usando i parametri GIA calibrati in ex2
% (calibrate_surface) e i dati di mercato gia bootstrappati. Non ripete piu
% bootstrap / ATM vol / calibrazione: li riceve in input.
%
% INPUT:
%   params : struct dei parametri calibrati (output di calibrate_surface)
%            .MA.alpha, .MA.beta
%            .AB.k,     .AB.eta
%            .GL.alpha, .GL.beta
%   market : struct dei dati di mercato/supporto (output di calibrate_surface)
%            .strikes, .calls, .puts, .expiries
%            .discount_factor, .forward, .R2, .sigma_ATM, .yf
%
% Se chiamata senza argomenti, esegue prima calibrate_surface (modalita
% standalone). La forward-start window [t1 = yf(2), t2 = yf(4)] e condivisa;
% ogni modello e poi prezzato con il proprio engine:
%   MA      -> Functions/   (FA_simulation, pricing_fwd_start_analytic/_MC)
%   AB, GL  -> UNIFIED Simulation/ engine: model_marginal_cf dispatch ->
%              ccdf_increment_FFT -> tail_adjustment -> simulate_from_cdf ->
%              price_fwd_start_MC. The only model-specific inputs are the
%              marginal CF (cf_AB / cf_GL from Distributions/) and the I_0
%              normalisation.
%
% All three legs apply the Lemma 2 forward rescaling (Forward.pdf). We compare
% the forward-start price (analytic / FFT reference) against Monte Carlo: MA
% across many strikes; AB and GL at the ATM forward start (K2 = 1), where the
% payoff reduces to max(W, 0).
%
% OUTPUT:
%   LA_results : struct array con i prezzi forward-start AB/GL (analytic + MC).

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
addpath("Simulation/Simulation_MA/")

% Single RNG seed for reproducibility
rng(1234);

% Standalone mode: if no calibrated inputs are passed, calibrate first.
if nargin < 2
    [params, market] = calibrate_surface(struct('verbose', false));
end

fprintf('=========================================================================\n');
fprintf('        LINEAR ADDITIVE MODELS (MA / AB / GL) - FORWARD-START ENGINE      \n');
fprintf('=========================================================================\n\n');

%% =========================================================================
% STEP 1-2: UNPACK MARKET DATA AND CALIBRATED PARAMETERS
% =========================================================================
fprintf('STEP 1-2: Using pre-calibrated curve, ATM vol and model parameters...\n');

% Market / support data (from ex2 calibrate_surface)
strikes         = market.strikes;
calls           = market.calls;
puts            = market.puts;
expiries        = market.expiries;
discount_factor = market.discount_factor;
forward         = market.forward;
sigma_ATM       = market.sigma_ATM;
yf              = market.yf;
nT              = numel(expiries);

fprintf('  -> Market data and ATM volatility loaded for %d maturities.\n\n', nT);

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

% 1. Model Parameters (calibrated in Exercise 2 -> params.MA = [alpha; beta])
alpha_MA = params.MA(1);
beta_MA  = params.MA(2);
gamma_MA = (1/alpha_MA) - (1/beta_MA);
fprintf('  -> MA: alpha = %.6f, beta = %.6f (from ex2 calibration)\n', alpha_MA, beta_MA);

% 2. Relevant dates / vols on the forward-start window
yf_fwd        = [yf(iT1); yf(iT2)];
sigma_ATM_fwd = [sigma_ATM(iT1); sigma_ATM(iT2)];
df_fwd        = discount_factor(iT2);
forward_fwd   = forward(iT2);
fwd_factor_MA = discount_factor(iT1) / discount_factor(iT2);   % Lemma 2 (Forward.pdf)

% 3. Integrated Volatility computation
C_MA   = 1 / ((1/alpha_MA) + (1/beta_MA));
I0     = I0_MA(params.MA);
sigmat = (sigma_ATM_fwd / I0) .* sqrt(yf_fwd);

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
[price_MC, CI_MC] = pricing_fwd_start_MC(forward_fwd, K2_vec, df_fwd, N_sim_MA, M_MA, dz_MA, sigmat, alpha_MA, beta_MA, fwd_factor_MA);

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

%% =========================================================================
% STEP 7 (MA): MOMENTS CHECK
% =========================================================================
compare_moments_MA(N_sim_MA, M_MA, dz_MA, sigmat, params.MA)

%% #########################################################################
%  MODELS 2 & 3 - ADDITIVE BACHELIER (AB) + GENERALIZED LOGISTIC (GL)
%  Unified Simulation/ engine. Both models share the EXACT same pipeline:
%    cdf (ccdf_increment_FFT) -> clean (tail_adjustment) ->
%    simulate (simulate_from_cdf) -> plot (plot_mc_check) ->
%    price (price_fwd_start_MC).
%  The only model-specific input is the parameter struct + I_0 normalisation.
%  #########################################################################
fprintf('STEP 7 (AB / GL): Unified forward-start pricing...\n');

N_sim_LA   = 1e6;
fwd_factor = discount_factor(iT1) / discount_factor(iT2);   % Lemma 2 (Forward.pdf)

% model configuration (params + I_0 normalisation)
% AB: I_0 from I0_AB (ex2). GL: I_0 = sqrt(2*pi)*E[zeta_+] (Baviera-Massaria Eq.14)
% Parameters calibrated in Exercise 2 -> params.AB / params.GL
kAB = params.AB(1);  eta_AB = params.AB(2);
alpha_GL = params.GL(1);  beta_GL = params.GL(2);
fprintf('  -> AB: k = %.6f, eta = %.6f | GL: alpha = %.6f, beta = %.6f (from ex2)\n', ...
        kAB, eta_AB, alpha_GL, beta_GL);
integrand_mean_GL = @(x) pdf_GL([alpha_GL; beta_GL], x) .* x;

% params is a column vector per the project-wide convention:
%   AB -> [k; eta] , GL -> [alpha; beta]   (see model_marginal_cf)
models(1) = struct('name', 'AB', ...
                   'params', [kAB; eta_AB], ...
                   'I0', I0_AB(0, [kAB; eta_AB]));
models(2) = struct('name', 'GL', ...
                   'params', [alpha_GL; beta_GL], ...
                   'I0', sqrt(2*pi) * quadgk(integrand_mean_GL, 0, inf));

LA_results = struct('name', {}, 'price_an', {}, 'price_mc', {}, 'IC', {});

for m = 1:numel(models)
    name   = models(m).name;
    params = models(m).params;
    
    sigma_t  = sigma_ATM / models(m).I0;
    sigma_T1 = sigma_t(iT1);
    sigma_T2 = sigma_t(iT2);

    % Evaluation grid scaled to the increment std (model-agnostic proxy):
    %   Var(W) ~ sigma_T2^2*T2 - fwd_factor^2*sigma_T1^2*T1.
    std_W  = sqrt(sigma_T2^2*T2 - fwd_factor^2*sigma_T1^2*T1);
    x_grid = linspace(-12*std_W, 12*std_W, 2000)';

    % conditional CDF t1 -> t2 (Lewis two-shift + Lemma 2) 
    cdf_raw        = ccdf_increment_FFT(name, params, T1, T2, sigma_T1, sigma_T2, x_grid, fwd_factor);
    [cdf_c, x_c]   = tail_adjustment(x_grid, cdf_raw, 10);

    % simulate the increment + visual check 
    Z_inc = simulate_from_cdf(cdf_c, x_c, true, N_sim_LA);
    plot_mc_check(Z_inc, x_c, cdf_c, T1, T2);
    set(gcf, 'Name', sprintf('%s increment t1 -> t2', name));

    % forward-start K2 = 1: MC (full path, Lemma 2)
    [price_mc, IC] = price_fwd_start_MC(name, params, T1, T2, sigma_T1, sigma_T2, ...
                        N_sim_LA, forward(iT2), discount_factor(iT1), ...
                        discount_factor(iT2), x_grid, 1);

    % FFT reference: E[max(W,0)] = int_0^inf (1 - F_W) dx
    F0       = interp1(x_c, cdf_c, 0, 'spline');
    mask_pos = x_c > 0;
    price_an = discount_factor(iT2) * trapz([0; x_c(mask_pos)], 1 - [F0; cdf_c(mask_pos)]);

    LA_results(m) = struct('name', name, 'price_an', price_an, ...
                           'price_mc', price_mc, 'IC', IC);

    fprintf('\n--- %s FORWARD START (K2 = 1) ---\n', name);
    fprintf('  FFT reference (analytic) : %.6f\n', price_an);
    fprintf('  Monte Carlo              : %.6f   95%% CI [%.6f, %.6f]\n', ...
            price_mc, IC(1), IC(2));
    fprintf('  Difference               : %.2f bps\n\n', abs(price_mc - price_an)*1e4);
end

%% =========================================================================
% STEP 8: UNIFORM FORWARD-START via lewis_FFT_call (all 3 models, K2 = 1)
% =========================================================================
fprintf('STEP 8: Uniform forward-start via lewis_FFT_call (AB/GL)...\n');
M_lewis  = 16;
dz_lewis = 0.05;
fwd_LA   = discount_factor(iT1) / discount_factor(iT2);   % Lemma 2 forward factor
% Use the per-model parameter scalars already extracted above (the loop above
% reassigns the local `params` to each model's vector, so params.* is gone).
par_AB = [kAB; eta_AB];  par_GL = [alpha_GL; beta_GL];
uspec = {'AB', par_AB, @(u,p,sc) cf_increment_AB(u,p,sc,fwd_LA), I0_AB(0, par_AB); ...
         'GL', par_GL, @(u,p,sc) cf_increment_GL(u,p,sc,fwd_LA), I0_GL(par_GL)};
price_uniform = zeros(size(uspec, 1), 1);
for i = 1:size(uspec, 1)
    nm = uspec{i,1}; par = uspec{i,2}; cfi = uspec{i,3}; I0u = uspec{i,4};
    scu = [(sigma_ATM(iT1)/I0u) * sqrt(T1), (sigma_ATM(iT2)/I0u) * sqrt(T2)];
    % z = 0 (ATM increment call); doubleshift = true for tail accuracy.
    price_uniform(i) = discount_factor(iT2) * ...
        lewis_FFT_call(cfi, M_lewis, dz_lewis, par, scu, 0, true, nm);
    fprintf('  -> %s : %.6f\n', nm, price_uniform(i));
end
fprintf('\n');

%% =========================================================================
% FINAL SUMMARY (K2 = 1)
% =========================================================================
fprintf('=========================================================================\n');
fprintf('              FORWARD-START SUMMARY AT K2 = 1 (ATM)                       \n');
fprintf('=========================================================================\n');
fprintf('%-6s | %-12s | %-12s | %-14s\n', 'Model', 'analytic', 'MC', 'uniform(Lewis)');
fprintf('-------------------------------------------------------------------------\n');
for m = 1:numel(LA_results)
    iu = find(strcmp(uspec(:,1), LA_results(m).name));
    fprintf('%-6s | %-12.6f | %-12.6f | %-14.6f\n', ...
            LA_results(m).name, LA_results(m).price_an, LA_results(m).price_mc, ...
            price_uniform(iu));
end
fprintf('=========================================================================\n');
fprintf('                 SIMULATION COMPLETED SUCCESSFULLY.                       \n');
fprintf('=========================================================================\n');
