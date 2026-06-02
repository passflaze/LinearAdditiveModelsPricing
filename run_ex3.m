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
% ogni modello e poi prezzato con il singolo engine UNIFICATO
% pricing_fwd_start_MC (switch su 'MA' | 'AB' | 'GL'):
%   MA      -> delega a pricing_fwd_start_MA_MC (FA_simulation: IA base +
%              FA increment).
%   AB, GL  -> CDF di entrambe le leg via lewis_FFT_digital sulla CF increment
%              Lemma-2 (cf_increment_AB / cf_increment_GL) -> simulate_from_cdf.
%              Gli unici input model-specific sono i parametri + la I_0.
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
addpath(genpath("Pricing/"));   % Pricing/Analytic (fwd-start AB/GL/MA) + Pricing/MC

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
fwd_factor    = discount_factor(iT1) / discount_factor(iT2);  % Lemma 2 (Forward.pdf) forward rescaling

% 3. Integrated Volatility computation
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

% Analytical Pricing (vectorized). fwd_factor propagates Lemma 2, so the
% analytic prices the SAME forward convention as the MC (consistent comparison).
price_analytic = pricing_fwd_start_analytic(alpha_MA, beta_MA, sigmat, df_fwd, K2_vec, forward_fwd, fwd_factor);

% Monte Carlo Pricing (vectorized) via the unified forward-start engine.
% sigmat already carries sqrt(T) (= the full scale factors at t1, t2).
rng(1234);
fprintf('  -> Running Monte Carlo Simulation (%d paths)...\n', N_sim_MA);
[price_MC, CI_MC] = pricing_fwd_start_MC('MA', params.MA, sigmat(1), sigmat(2), ...
    0, M_MA, dz_MA, forward_fwd, discount_factor(iT1), discount_factor(iT2), K2_vec);

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

    % Lemma 2: the increment s leg is rescaled by fwd_factor (same convention as
    % the MC / analytic) so this exact-CDF cross-check is consistent at ATM.
    sigmat_s = fwd_factor * sigmat(1);
    ps_plus  = beta_MA / sigmat_s;
    ps_minus = alpha_MA / sigmat_s;
    pt_plus  = beta_MA / sigmat(2);
    pt_minus = alpha_MA / sigmat(2);
    drift_t1_t2 = gamma_MA * (sigmat(2) - sigmat_s);

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
%  Unified engine pricing_fwd_start_MC (case 'AB'/'GL'). Both models share the
%  EXACT same pipeline inside it:
%    cdf (lewis_FFT_digital) -> clean (tail_adjustment, internal) ->
%    simulate (simulate_from_cdf) -> price.
%  The only model-specific input is the parameter struct + I_0 normalisation.
%  #########################################################################
fprintf('STEP 7 (AB / GL): Unified forward-start pricing...\n');

N_sim_LA   = 1e6;
M_lewis    = 16;        % Lewis-FFT grid exponent (N = 2^M) for AB/GL CDFs
dz_lewis   = 0.05;      % Lewis-FFT z-step (dollar increments)

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
    
    % FULL scale factors at t1, t2 (= (sigma_ATM/I0)*sqrt(T)).
    sigma_t  = sigma_ATM / models(m).I0;
    sc_T1    = sigma_t(iT1) * sqrt(T1);
    sc_T2    = sigma_t(iT2) * sqrt(T2);

    % forward-start K2 = 1: MC (full path, Lemma 2). The unified engine builds
    % both legs' CDFs via lewis_FFT_digital and returns the conditional CDF
    % (diag) for the visual check / FFT reference.
    rng(1234);
    [price_mc, IC, diag] = pricing_fwd_start_MC(name, params, sc_T1, sc_T2, ...
                        0, M_lewis, dz_lewis, forward(iT2), ...
                        discount_factor(iT1), discount_factor(iT2), 1);

    % visual check on the conditional increment t1 -> t2
    x_c = diag.x_cond;  cdf_c = diag.cdf_cond;
    plot_mc_check(diag.W, x_c, cdf_c, T1, T2);
    set(gcf, 'Name', sprintf('%s increment t1 -> t2', name));

    % Semi-analytic forward-start (K2 = 1) via Lewis-FFT on the Lemma-2
    % increment (Pricing/Analytic). At K2 = 1 this is df * E[max(W,0)].
    if strcmp(name, 'AB')
        price_an = pricing_fwd_start_AB_analytic(params, sc_T1, sc_T2, ...
                       discount_factor(iT2), 1, forward(iT2), fwd_factor, M_lewis, dz_lewis);
    else
        price_an = pricing_fwd_start_GL_analytic(params, sc_T1, sc_T2, ...
                       discount_factor(iT2), 1, forward(iT2), fwd_factor, M_lewis, dz_lewis);
    end

    LA_results(m) = struct('name', name, 'price_an', price_an, ...
                           'price_mc', price_mc, 'IC', IC);

    fprintf('\n--- %s FORWARD START (K2 = 1) ---\n', name);
    fprintf('  Analytic (Lewis-FFT)     : %.6f\n', price_an);
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
