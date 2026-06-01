function compare_moments_MA(N_sim, M, dz, scale_factor, params)
% COMPARE_MOMENTS_MA Computes and compares empirical (MC) vs analytical 
% moments for the Minimal Additive (MA) finite activity increment.
%
% Inputs:
%   forward, K, df : Standard pricing parameters (passed for signature consistency)
%   N_sim          : Number of Monte Carlo paths
%   M, dz          : Grid parameters for FFT/simulation
%   scale_factor   : Integrated volatilities [sigma*sqrt(t1), sigma*sqrt(t2)]
%   params         : [alpha_MA, beta_MA]

    % =========================================================================
    % STEP 1: INTERNAL PARAMETER COMPUTATION
    % =========================================================================
    alpha_MA = params(1);
    beta_MA = params(2);
    gamma_MA = (1/alpha_MA) - (1/beta_MA);
    
    % Tail decay parameters using scale_factor
    ps_plus  = beta_MA  / scale_factor(1);
    ps_minus = alpha_MA / scale_factor(1);
    
    pt_plus  = beta_MA  / scale_factor(2);
    pt_minus = alpha_MA / scale_factor(2);
    
    % Drift vector calculation (from t1 to t2)
    drift_t1_t2 = gamma_MA * (scale_factor(2) - scale_factor(1));

    % =========================================================================
    % STEP 2: MONTE CARLO SIMULATION & EMPIRICAL MOMENTS
    % =========================================================================
    fprintf('\n========================================================\n');
    fprintf('  MA MOMENTS COMPARISON: MONTE CARLO vs ANALYTICAL\n');
    fprintf('========================================================\n');
    
    fprintf('Running Monte Carlo simulation with %d paths...\n', N_sim);
    tic;
    % Finite-activity increment [t1 -> t2]: the 2-vector scale_factor
    % [sigma_t1; sigma_t2] is required so conditional_cf_MA_FA gets both tails.
    % params = [alpha; beta]; delta_mu = drift_t1_t2 is added in the finite branch.
    ft1t2 = FA_simulation(N_sim, M, dz, drift_t1_t2, ...
                          pt_plus, pt_minus, ps_plus, ps_minus, 1, 'finite', 1, ...
                          params(:), scale_factor(:), []);
    time_MC = toc;
    
    % Empirical Moments (Corrected: Mean, Variance, Skewness, Kurtosis)
    M1_MC = mean(ft1t2);
    M2_MC = var(ft1t2);      % Match analytical Variance, not STD
    M3_MC = skewness(ft1t2);
    M4_MC = kurtosis(ft1t2); % Pearson standard Kurtosis

    % =========================================================================
    % STEP 3: ANALYTICAL MOMENTS
    % =========================================================================
    tic;
    [M1_an, M2_an, M3_an, M4_an] = moments_generator(pt_plus, pt_minus, ...
                                                     ps_plus, ps_minus, drift_t1_t2);
    time_an = toc;

    % =========================================================================
    % STEP 4: DIAGNOSTIC PRINTING
    % =========================================================================
    % Absolute error for mean (since it clusters around 0), Relative for others
    err_M1 = abs(M1_MC - M1_an);
    err_M2 = abs(M2_MC - M2_an) / abs(M2_an) * 100;
    err_M3 = abs(M3_MC - M3_an) / abs(M3_an) * 100;
    err_M4 = abs(M4_MC - M4_an) / abs(M4_an) * 100;

    fprintf('\nCOMPUTATIONAL TIME:\n');
    fprintf('  Monte Carlo Time : %10.6f seconds\n', time_MC);
    fprintf('  Analytical Time  : %10.6f seconds\n', time_an);
    fprintf('  Speedup Factor   : %10.0fx faster\n', time_MC / time_an);
    
    fprintf('\nSTATISTICAL MOMENTS:\n');
    fprintf('  %-10s | %-12s | %-12s | %-12s\n', 'Moment', 'Analytical', 'Monte Carlo', 'Error');
    fprintf('  --------------------------------------------------------\n');
    fprintf('  %-10s | %12.6f | %12.6f | %10.2e (Abs)\n', '1. Mean', M1_an, M1_MC, err_M1);
    fprintf('  %-10s | %12.6f | %12.6f | %11.4f %%\n', '2. Variance', M2_an, M2_MC, err_M2);
    fprintf('  %-10s | %12.6f | %12.6f | %11.4f %%\n', '3. Skewness', M3_an, M3_MC, err_M3);
    fprintf('  %-10s | %12.6f | %12.6f | %11.4f %%\n', '4. Kurtosis', M4_an, M4_MC, err_M4);
    fprintf('========================================================\n\n');

    % =========================================================================
    % STEP 5: PLOTTING
    % =========================================================================
    figure('Name', 'MA Increment: Distribution and Moments', 'Position', [100, 100, 1000, 450]);
    
    % Subplot 1: Distribution Histogram
    subplot(1, 2, 1);
    histogram(ft1t2, 100, 'Normalization', 'pdf', 'FaceColor', [0.2 0.4 0.6], 'EdgeColor', 'none');
    hold on;
    xline(M1_an, 'r--', 'LineWidth', 2, 'DisplayName', 'Analytical Mean');
    title('Simulated Distribution of \Delta f');
    xlabel('Value');
    ylabel('Probability Density');
    legend('Location', 'best');
    grid on;
    
    % Subplot 2: Relative Errors
    subplot(1, 2, 2);
    errors = [err_M2, err_M3, err_M4]; % Excluded mean to avoid scale distortion
    labels = {'Variance', 'Skewness', 'Kurtosis'};
    b = bar(errors, 'FaceColor', [0.8 0.3 0.3]);
    set(gca, 'xticklabel', labels);
    title('Monte Carlo Relative Error (%)');
    ylabel('Error (%)');
    ylim([0, max(errors) * 1.2 + 1e-4]); % Add some headroom
    
    % Add error values on top of bars
    for i = 1:length(errors)
        text(i, errors(i), sprintf('%.2f%%', errors(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    grid on;
end