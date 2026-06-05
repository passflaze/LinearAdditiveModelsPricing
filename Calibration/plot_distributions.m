function plot_distributions(eta_AB, kappa_AB, alpha_MA, beta_MA, alpha_GL, beta_GL)
% PLOT_DISTRIBUTIONS  Compare the MA, AB (NIG) and GL densities vs a standard
%   normal, on linear and semi-log scales.
%
% INPUTS:
%   eta_AB, kappa_AB   : AB skew and IG parameters
%   alpha_MA, beta_MA  : MA left/right rates
%   alpha_GL, beta_GL  : GL left/right shape parameters
% OUTPUT:
%   none (draws a linear-scale and a semi-log-scale density figure)

    x_grid = linspace(-30, 30, 500);

    % MA density (sigmat = 1, no time-scaling).
    y_MA = pdf_MA([alpha_MA; beta_MA], 1, x_grid);

    % AB density via the NIG parametrization (x, alpha, beta, mu, delta).
    mu        = eta_AB;
    gamma_AB  = -eta_AB;
    lam_IG1   = 1 / kappa_AB;
    alpha_nig = sqrt(gamma_AB^2 + lam_IG1);
    beta_nig  = gamma_AB;
    delta_nig = sqrt(lam_IG1);
    mu_nig    = mu;
    y_AB = pdf_NIG(x_grid, alpha_nig, beta_nig, mu_nig, delta_nig);

    % GL density.
    y_GL = pdf_GL([alpha_GL; beta_GL], x_grid);

    % Standard normal benchmark.
    mu_norm = 0;
    sigma_norm = 1;
    y_norm = (1 / (sigma_norm * sqrt(2*pi))) .* exp(-0.5 .* ((x_grid - mu_norm) ./ sigma_norm).^2);

    % --- Plot 1: linear scale ---
    figure('Color', 'w', 'Position', [100, 100, 900, 600]);
    hold on; grid on; box on;

    plot(x_grid, y_MA, 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Asymmetric Laplace (MA)');
    plot(x_grid, y_AB, 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB Distribution');
    plot(x_grid, y_GL, 'LineWidth', 2.5, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'Generalized Logistic (GL)');
    plot(x_grid, y_norm, '--k', 'LineWidth', 2, 'DisplayName', 'Standard Normal Benchmark');

    title('Probability Density Functions (Linear Scale)', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Density', 'FontSize', 14);
    legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');
    xlim([min(x_grid), max(x_grid)]);
    ylim([0, max([y_MA, y_AB, y_GL, y_norm]) * 1.1]);
    hold off;

    % --- Plot 2: semi-log scale ---
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

