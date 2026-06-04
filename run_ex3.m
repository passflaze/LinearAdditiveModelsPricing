function LA_results = run_ex3(params, market, opts)
%RUN_EX3 Computes Forward Start prices using Linear Additive Models 
%        (GL, MA, AB) and compares Analytic vs Monte Carlo results.

    % =========================================================================
    % 1. INITIALIZATION & MARKET SETUP
    % =========================================================================
    addpath("Utilities/");
    addpath("Distributions/");
    addpath("Calibration/");
    addpath("Calibration/Calibration_AB/");
    addpath("Calibration/Calibration_MA/");
    addpath("Calibration/Calibration_GL/");
    addpath("Simulation/");
    addpath("Simulation/Simulation_MA/");
    addpath(genpath("Pricing/"));
    
    % --- Options Initialization ---
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'verbose'), opts.verbose = false; end
    if ~isfield(opts, 'plot'),    opts.plot = true;     end
    % ------------------------------

    % Standalone mode: if no calibrated inputs are passed, calibrate first
    if nargin < 2
        [params, market] = calibrate_surface(struct('verbose', false));
    end
    
    fprintf('=========================================================================\n');
    fprintf('        LINEAR ADDITIVE MODELS (MA / AB / GL) - FORWARD-START ENGINE     \n');
    fprintf('=========================================================================\n\n');
    fprintf('STEP 1: Using pre-calibrated curve, ATM vol and model parameters...\n');
    
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

    % =========================================================================
    % 2. GLOBAL PRICING PARAMETERS
    % =========================================================================
    iT1 = 2;  
    iT2 = 4;
    T1  = yf(iT1);  
    T2  = yf(iT2);
    K2 = opts.K2;
    fprintf('  Forward-start window: t1 = %.4f y (idx %d), t2 = %.4f y (idx %d)\n\n', T1, iT1, T2, iT2);
    
    fwd_factor = discount_factor(iT1) / discount_factor(iT2);

    % =========================================================================
    % 3. FFT VALIDATION (GAUSSIAN CDF SANITY CHECK)
    % =========================================================================
    if opts.verbose
        fprintf('--- 0. FFT VALIDATION ---\n');
        M_gauss     = 14;
        dz_gauss    = 2.5e-3;
        shift_gauss = 0.5;
        [cdf_clean_gauss, z_grid_gauss] = lewis_fft_algorithm_trial(M_gauss, dz_gauss, shift_gauss, 1);
        cdf_analytical = normcdf(z_grid_gauss);
        
        figure('Name', 'FFT vs Analytical Gaussian CDF', 'Color', 'white');
        plot(z_grid_gauss, cdf_clean_gauss, 'b-', 'LineWidth', 2, 'DisplayName', 'FFT Inversion');
        hold on;
        plot(z_grid_gauss, cdf_analytical, 'r--', 'LineWidth', 2, 'DisplayName', 'Analytical CDF');
        title('FFT Inversion vs Analytical Gaussian CDF');
        xlabel('Spatial Grid (z)'); ylabel('Cumulative Probability');
        legend('Location', 'northwest'); xlim([-5, 5]); grid on; hold off;
        
        fprintf('  -> FFT inversion validation against Gaussian CDF complete.\n\n');
    end

    % =========================================================================
    % 4. MODEL PRICING EXECUTION
    % =========================================================================

    % --- 4.1 GENERALIZED LOGISTIC (GL) MODEL ---
    fprintf('--- 1. GENERALIZED LOGISTIC (GL) MODEL ---\n');
    params_GL = params.GL;
    I0_GL_val = I0_GL(params_GL);
    
    sigma_t_GL = sigma_ATM / I0_GL_val;
    sc_T1_GL   = sigma_t_GL(iT1) * sqrt(T1);
    sc_T2_GL   = sigma_t_GL(iT2) * sqrt(T2);

    M_GL = opts.mc.M_GL;
    dz_GL = opts.mc.dz_GL;
    Nsim_GL = opts.mc.Nsim_GL;
    
    rng(1234);
    fprintf('  > Computing Forward Start ATM (MC)...\n');
    [price_mc_GL, IC_GL, diag_GL] = pricing_fwd_start_MC('GL', params_GL, sc_T1_GL, sc_T2_GL, ...
                        Nsim_GL, M_GL, dz_GL, forward(iT2), discount_factor(iT1), discount_factor(iT2), K2, opts);
                        
    if opts.plot
        plot_mc_check(diag_GL.W, diag_GL.x_cond, diag_GL.cdf_cond, T1, T2);
        set(gcf, 'Name', 'GL increment t1 -> t2', 'Color', 'white');
    end
    if K2 ==1
    fprintf('  > Computing Forward Start ATM (Analytic via Lewis-FFT)...\n');
    price_an_GL = pricing_fwd_start_GL_analytic(params_GL, sc_T1_GL, sc_T2_GL, ...
                       discount_factor(iT2), K2, forward(iT2), fwd_factor, M_GL, dz_GL);
                       
    fprintf('  > Computing Uniform Forward-Start via lewis_FFT_call...\n');
    cfi_GL = @(u,p,sc) cf_increment_GL(u,p,sc,fwd_factor);
    price_uniform_GL = discount_factor(iT2) * lewis_FFT_call(cfi_GL, M_GL, dz_GL, params_GL, [sc_T1_GL, sc_T2_GL], 0, true, 'GL');

    fprintf('  [Results GL] Analytic: %.6f | MC: %.6f | Uniform: %.6f\n\n', price_an_GL, price_mc_GL, price_uniform_GL);
    end

    % --- 4.2 MINIMAL ADDITIVE (MA) MODEL ---
    fprintf('--- 2. MINIMAL ADDITIVE (MA) MODEL ---\n');
    params_MA = params.MA;
    alpha_MA  = params_MA(1);
    beta_MA   = params_MA(2);
    gamma_MA  = (1/alpha_MA) - (1/beta_MA);
    
    I0_MA_val = I0_MA(params_MA);
    sigma_t_MA = sigma_ATM / I0_MA_val;
    sc_T1_MA   = sigma_t_MA(iT1) * sqrt(T1);
    sc_T2_MA   = sigma_t_MA(iT2) * sqrt(T2);
    sigmat_MA  = [sc_T1_MA; sc_T2_MA];

    M_MA = opts.mc.M_MA;
    dz_MA = opts.mc.dz_MA;
    Nsim_MA = opts.mc.Nsim_MA;
    
    % Multi-Strike Setup
    num_random_strikes = 20;
    K2_vec = sort([K2, 0.8 + 0.4 * rand(1, num_random_strikes)]);
    idx_atm_MA = find(K2_vec == K2, 1);
    
    fprintf('  > Computing Forward Start Multi-Strike (Analytic)...\n');
    price_an_MA_vec = pricing_fwd_start_MA_analytic(alpha_MA, beta_MA, sigmat_MA, discount_factor(iT2), K2_vec, forward(iT2), fwd_factor);
    
    rng(1234);
    fprintf('  > Computing Forward Start Multi-Strike (MC)...\n');
    [price_mc_MA_vec, IC_MA_vec] = pricing_fwd_start_MC('MA', params_MA, sc_T1_MA, sc_T2_MA, ...
        Nsim_MA, M_MA, dz_MA, forward(iT2), discount_factor(iT1), discount_factor(iT2), K2_vec, opts);
        
    diff_bps_MA = (price_mc_MA_vec - price_an_MA_vec) * 10000;
    
    if opts.plot
        figure('Name', 'MA Forward Start Pricing', 'Color', 'white');
        hold on; grid on;
        plot(K2_vec, price_an_MA_vec, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytic Price');
        plot(K2_vec, price_mc_MA_vec, 'ro', 'MarkerSize', 6, 'DisplayName', 'Monte Carlo Price');
        if ~isempty(idx_atm_MA)
            plot(K2_vec(idx_atm_MA), price_an_MA_vec(idx_atm_MA), 'kp', 'MarkerSize', 15, 'MarkerFaceColor', 'y', 'DisplayName', 'ATM Forward (K_2 = 1)');
        end
        title('MA Forward Start Pricing: Analytic vs Monte Carlo');
        xlabel('Proportional Strike (K_2)'); ylabel('Discounted Option Price');
        legend('Location', 'best'); xlim([min(K2_vec), max(K2_vec)]); hold off;
    end
    
    % Extract ATM values for final summary table
    if ~isempty(idx_atm_MA)
        price_an_MA = price_an_MA_vec(idx_atm_MA);
        price_mc_MA = price_mc_MA_vec(idx_atm_MA);
        
        fprintf('  > Verifying ATM Analytical Price via Numerical Integration...\n');
        sigmat_s = fwd_factor * sc_T1_MA;
        drift_t1_t2 = gamma_MA * (sc_T2_MA - sigmat_s);
        [cdf_exact, x_grid] = exact_ma_increment_cdf(beta_MA/sc_T2_MA, alpha_MA/sc_T2_MA, beta_MA/sigmat_s, alpha_MA/sigmat_s, drift_t1_t2, opts.plot);
        cdf_fun = @(x) interp1(x_grid, cdf_exact, x, 'pchip', 'extrap');
        price_integrated_MA = discount_factor(iT2) * integral(@(x) 1 - cdf_fun(x), 0, x_grid(end));
        
        fprintf('    [Results MA] Formula: %.6f | Integral: %.6f | Error: %.2f bps\n', price_an_MA, price_integrated_MA, abs(price_integrated_MA - price_an_MA)*10000);
    else
        price_an_MA = NaN; price_mc_MA = NaN;
    end
    
    fprintf('  > Checking Moments...\n');
    compare_moments_MA(1e7, M_MA, dz_MA, sigmat_MA, params_MA, opts.plot);
    fprintf('\n');

    % --- 4.3 ADDITIVE BACHELIER (AB) MODEL ---
    fprintf('--- 3. ADDITIVE BACHELIER (AB) MODEL ---\n');
    params_AB = params.AB;
    I0_AB_val = I0_AB(0, params_AB);
    
    sigma_t_AB = sigma_ATM / I0_AB_val;
    sc_T1_AB   = sigma_t_AB(iT1) * sqrt(T1);
    sc_T2_AB   = sigma_t_AB(iT2) * sqrt(T2);
    dz_AB = opts.mc.dz_AB;
    M_AB = opts.mc.M_AB;
    Nsim_AB = opts.mc.Nsim_AB;
    
    rng(1234);
    fprintf('  > Computing Forward Start ATM (MC)...\n');
    [price_mc_AB, IC_AB, diag_AB] = pricing_fwd_start_MC('AB', params_AB, sc_T1_AB, sc_T2_AB, ...
                        Nsim_AB, M_AB, dz_AB, forward(iT2), discount_factor(iT1), discount_factor(iT2), K2, opts);
                        
    if opts.plot
        plot_mc_check(diag_AB.W, diag_AB.x_cond, diag_AB.cdf_cond, T1, T2);
        set(gcf, 'Name', 'AB increment t1 -> t2', 'Color', 'white');
    end
    if K2==1
        fprintf('  > Computing Forward Start ATM (Analytic via Lewis-FFT)...\n');
        price_an_AB = pricing_fwd_start_AB_analytic(params_AB, sc_T1_AB, sc_T2_AB, ...
                           discount_factor(iT2), K2, forward(iT2), fwd_factor, M_AB, dz_AB);
                           
        fprintf('  > Computing Uniform Forward-Start via lewis_FFT_call...\n');
        cfi_AB = @(u,p,sc) cf_increment_AB(u,p,sc,fwd_factor);
        price_uniform_AB = discount_factor(iT2) * lewis_FFT_call(cfi_AB, M_AB, dz_AB, params_AB, [sc_T1_AB, sc_T2_AB], 0, true, 'AB');
    
        fprintf('  [Results AB] Analytic: %.6f | MC: %.6f | Uniform: %.6f\n\n', price_an_AB, price_mc_AB, price_uniform_AB);
    end

   % =========================================================================
    % 5. RESULTS AGGREGATION & REPORTING
    % =========================================================================
    models_names = {'GL'; 'MA'; 'AB'};
    
    if K2 == 1
        prices_an = [price_an_GL; price_an_MA; price_an_AB];
        prices_mc = [price_mc_GL; price_mc_MA; price_mc_AB];
        prices_un = [price_uniform_GL; NaN; price_uniform_AB]; % MA does not use uniform lewis_FFT_call
        
        % Calculate the difference in basis points
        diff_bps = (prices_mc - prices_an) * 10000;
        
        SummaryTable = table(models_names, prices_an, prices_mc, diff_bps, prices_un, ...
            'VariableNames', {'Model', 'Analytic', 'MonteCarlo', 'Diff_bps', 'Uniform_Lewis'});
            
        fprintf('=========================================================================\n');
        fprintf('              FORWARD-START SUMMARY AT K2 = 1 (ATM)                      \n');
        fprintf('=========================================================================\n');
        disp(SummaryTable);
    else
        prices_mc = [price_mc_GL; price_mc_MA; price_mc_AB];
        
        % Extract Confidence Intervals for the specific K2 strike
        IC_MA = IC_MA_vec(:, idx_atm_MA);
        
        CI_Lower = [IC_GL(1); IC_MA(1); IC_AB(1)];
        CI_Upper = [IC_GL(2); IC_MA(2); IC_AB(2)];
        
        SummaryTable = table(models_names, prices_mc, CI_Lower, CI_Upper, ...
            'VariableNames', {'Model', 'MonteCarlo', 'CI_Lower', 'CI_Upper'});
            
        fprintf('=========================================================================\n');
        fprintf('              FORWARD-START SUMMARY AT K2 = %.4f                         \n', K2);
        fprintf('=========================================================================\n');
        disp(SummaryTable);
    end
    
    fprintf('=========================================================================\n');
    fprintf('                 SIMULATION COMPLETED SUCCESSFULLY.                      \n');
    fprintf('=========================================================================\n');
    
    % Output Struct
    LA_results.PricingData = SummaryTable;
end