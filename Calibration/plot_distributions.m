function plot_distributions(eta_AB, kappa_AB, alpha_MA, beta_MA, alpha_GL, beta_GL)
% COMPARE_DISTRIBUTIONS Analyzes and compares the analytical PDFs of three 
% quantitative finance distributions (MA, AB, GL) against a Standard Normal.
%
% INPUTS:
%   eta_AB   : Drift/skew parameter for the Additive Bachelier (NIG) model
%   kappa_AB : Kurtosis/variance parameter for the Additive Bachelier (NIG) model
%   alpha_MA : Left tail shape parameter for the Asymmetric Laplace model
%   beta_MA  : Right tail shape parameter for the Asymmetric Laplace model
%   alpha_GL : Left tail shape parameter for the Generalized Logistic model
%   beta_GL  : Right tail shape parameter for the Generalized Logistic model

    % =========================================================================
    % 1. Global Parameters and Evaluation Grid
    % =========================================================================
    % Unified spatial grid for evaluating the probability density functions
    x_grid = linspace(-30, 30, 500);

    % =========================================================================
    % 2. Asymmetric Laplace Distribution (MA Model)
    % =========================================================================
    % Evaluate MA density using dedicated function from Distributions folder
    % Note: sigmat=1 (no time-scaling in this visualization context)
    y_MA = pdf_MA(alpha_MA, beta_MA, 1, x_grid);

    % =========================================================================
    % 3. Additive Bachelier Distribution (AB via NIG Model)
    % =========================================================================
    % Structural parameter mapping for the Normal Inverse Gaussian
    mu        = eta_AB;
    gamma_AB  = -eta_AB;
    lam_IG1   = 1 / kappa_AB;

    alpha_nig = sqrt(gamma_AB^2 + lam_IG1);
    beta_nig  = gamma_AB;
    delta_nig = sqrt(lam_IG1);
    mu_nig    = mu;

    % Evaluate AB density using custom NIG PDF function from Distributions folder
    % Note: parameter order is (x, alpha, beta, mu, delta) vs. MATLAB's nigpdf
    y_AB = pdf_NIG(x_grid, alpha_nig, beta_nig, mu_nig, delta_nig);

    % =========================================================================
    % 4. Generalized Logistic Distribution (GL)
    % =========================================================================
    % Evaluate GL density using dedicated function from Distributions folder
    y_GL = pdf_GL(alpha_GL, beta_GL, x_grid);

    % =========================================================================
    % 5. Benchmark: Standard Normal Distribution
    % =========================================================================
    mu_norm = 0;
    sigma_norm = 1;
    % Explicit formula to avoid Statistics Toolbox dependencies
    y_norm = (1 / (sigma_norm * sqrt(2*pi))) .* exp(-0.5 .* ((x_grid - mu_norm) ./ sigma_norm).^2);

    % =========================================================================
    % PLOT 1: Linear Scale Comparison (vs. Standard Normal)
    % =========================================================================
    figure('Color', 'w', 'Position', [100, 100, 900, 600]);
    hold on; grid on; box on;
    
    % Plotting main distributions with professional MATLAB color palette
    plot(x_grid, y_MA, 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Asymmetric Laplace (MA)');
    plot(x_grid, y_AB, 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB Distribution');
    plot(x_grid, y_GL, 'LineWidth', 2.5, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'Generalized Logistic (GL)');
    
    % Adding the Normal curve as a dashed benchmark
    plot(x_grid, y_norm, '--k', 'LineWidth', 2, 'DisplayName', 'Standard Normal Benchmark');
    
    % Formatting and aesthetics
    title('Probability Density Functions (Linear Scale)', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Density', 'FontSize', 14);
    legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');
    
    % Axis limits
    xlim([min(x_grid), max(x_grid)]);
    ylim([0, max([y_MA, y_AB, y_GL, y_norm]) * 1.1]);
    hold off;

    % =========================================================================
    % PLOT 2: Global Semi-Logarithmic Scale
    % =========================================================================
    figure('Color', 'w', 'Position', [150, 120, 900, 600]);
    hold on; grid on; box on;
    
    semilogy(x_grid, y_MA, 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Asymmetric Laplace (MA)');
    semilogy(x_grid, y_AB, 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB Distribution');
    semilogy(x_grid, y_GL, 'LineWidth', 2.5, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'Generalized Logistic (GL)');
    
    set(gca, 'YScale', 'log');
    
    title('Probability Density Functions (Semi-Log Scale)', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Density (Log)', 'FontSize', 14);
    legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');
    
    xlim([min(x_grid), max(x_grid)]);
    ylim([1e-10, max([y_MA, y_AB, y_GL]) * 1.1]);
    hold off;

    
    
end