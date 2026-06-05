function run_ex0(check_plot)
% RUN_EX0  Analytical PDF comparison of MA, Additive Bachelier, and GL models.
%   Draws linear-scale, semi-log, and tail-zoom PDF plots.
%   Optional Monte Carlo check overlays the empirical AB density.
%
% INPUTS:
%   check_plot - (logical, default false) run 1M AB Monte Carlo check

    if nargin < 1
        check_plot = false;
    end

    addpath('Distributions');

    %% Global parameters
    p_plus  = 1.5;
    p_minus = 0.9;

    x_grid = linspace(-30, 30, 2000);

    %% MA PDF (Asymmetric Laplace)
    % sigma*sqrt(t) = 1 for unit-time evaluation
    y_MA = pdf_MA([p_minus; p_plus], 1, x_grid);

    %% AB PDF (NIG)
    kappa_AB = 1 / (p_plus * p_minus);
    eta_AB   = ((1/kappa_AB) - (p_minus^2)) / (2 * p_minus);

    % NIG parameter mapping
    mu_nig    = eta_AB;
    gamma_AB  = -eta_AB;
    mu_IG1    = 1;
    lam_IG1   = 1 / kappa_AB;
    alpha_nig = sqrt(gamma_AB^2 + lam_IG1);
    beta_nig  = gamma_AB;
    delta_nig = sqrt(lam_IG1);

    y_AB = pdf_NIG(x_grid, alpha_nig, beta_nig, mu_nig, delta_nig);

    %% GL PDF
    y_GL = pdf_GL([p_minus; p_plus], x_grid);

    %% AB Monte Carlo check (optional)
    if check_plot
        numSim = 1e6;
        W        = random('InverseGaussian', mu_IG1, lam_IG1, [numSim, 1]);
        Z        = randn(numSim, 1);
        y_sim_AB = mu_nig + W .* gamma_AB + sqrt(W) .* Z;

        [f_emp_AB, ~] = ksdensity(y_sim_AB, x_grid);

        figure('Color', 'w', 'Position', [150, 150, 900, 600]);
        hold on;
        grid on;
        box on;
        plot(x_grid, f_emp_AB, 'LineWidth', 4, 'Color', [0.4940 0.1840 0.5560], ...
             'DisplayName', 'Empirical PDF (Monte Carlo 1M)');
        plot(x_grid, y_AB, '--', 'LineWidth', 2.5, 'Color', [0.9290 0.6940 0.1250], ...
             'DisplayName', 'Analytical PDF (NIG)');
        
        title('Additive Bachelier: Monte Carlo Simulation vs Analytical PDF', 'FontSize', 15, 'FontWeight', 'bold');
        xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('Density', 'FontSize', 14);
        legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');
        xlim([-8, 8]);
        ylim([0, max(y_AB) * 1.1]);
        hold off;
    end

    %% PDF comparison (linear scale)
    y_norm = (1 / sqrt(2*pi)) .* exp(-0.5 .* x_grid.^2);

    figure('Color', 'w', 'Position', [100, 100, 900, 600]);
    hold on;
    grid on;
    box on;
    plot(x_grid, y_MA, 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410],      'DisplayName', 'Asymmetric Laplace (MA)');
    plot(x_grid, y_AB, 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980],  'DisplayName', 'AB Distribution');
    plot(x_grid, y_GL, 'LineWidth', 2.5, 'Color', [0.4940 0.1840 0.5560],  'DisplayName', 'Generalized Logistic (GL)');
    plot(x_grid, y_norm, '--k', 'LineWidth', 2,                             'DisplayName', 'Normal Benchmark');
    
    title('Probability Density Functions Comparison (Linear Scale)', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Density', 'FontSize', 14);
    legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');
    
    xlim([-30, 30]);
    ylim([0, max([y_MA, y_AB, y_GL, y_norm]) * 1.1]);
    hold off;

    %% PDF comparison (semi-log scale)
    figure('Color', 'w', 'Position', [100, 100, 900, 600]);
    hold on;
    grid on;
    box on;
    semilogy(x_grid, y_MA, 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410],     'DisplayName', 'Asymmetric Laplace (MA)');
    semilogy(x_grid, y_AB, 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB Distribution');
    semilogy(x_grid, y_GL, 'LineWidth', 2.5, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'Generalized Logistic (GL)');
    set(gca, 'YScale', 'log');

    title('Probability Density Functions Comparison (Semi-Log Scale)', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Density', 'FontSize', 14);
    legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');

    xlim([-30, 30]);
    ylim([1e-10, max([y_MA, y_AB, y_GL]) * 1.1]);
    hold off;

    %% Tail analysis
    figure('Color', 'w', 'Position', [150, 150, 1000, 450]);

    x_left_tail  = linspace(-500, -450, 100);
    x_right_tail = linspace(450, 500, 100);

    y_GL_left  = pdf_GL([p_minus; p_plus], x_left_tail);
    y_MA_left  = pdf_MA([p_minus; p_plus], 1, x_left_tail);
    y_AB_left  = pdf_NIG(x_left_tail, alpha_nig, beta_nig, mu_nig, delta_nig);

    y_GL_right = pdf_GL([p_minus; p_plus], x_right_tail);
    y_MA_right = pdf_MA([p_minus; p_plus], 1, x_right_tail);
    y_AB_right = pdf_NIG(x_right_tail, alpha_nig, beta_nig, mu_nig, delta_nig);

    % Extreme Left Tail [-500, -450]
    subplot(1, 2, 1);
    hold on; grid on; box on;
    semilogy(x_left_tail, y_MA_left, 'LineWidth', 2, 'Color', [0 0.4470 0.7410],     'DisplayName', 'MA');
    semilogy(x_left_tail, y_AB_left, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB');
    semilogy(x_left_tail, y_GL_left, 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'GL');
    
    x_refL  = linspace(-500, -450, 50);
    anchorL = interp1(x_left_tail, y_AB_left, -450); % Anchored to the inner edge of the view
    refL    = anchorL * exp(p_minus * (x_refL + 450));
    semilogy(x_refL, refL, 'k--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('slope p_- = %.1f', p_minus));
        
    title('Extreme Left Tail (Log-Scale)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Log-Density', 'FontSize', 12);
    legend('Location', 'southwest', 'FontSize', 10);
    set(gca, 'YScale', 'log');
    xlim([-500, -450]);
    hold off;

    % Extreme Right Tail [450, 500]
    subplot(1, 2, 2);
    hold on; grid on; box on;
    semilogy(x_right_tail, y_MA_right, 'LineWidth', 2, 'Color', [0 0.4470 0.7410],     'DisplayName', 'MA');
    semilogy(x_right_tail, y_AB_right, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB');
    semilogy(x_right_tail, y_GL_right, 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'GL');
    
    x_refR  = linspace(450, 500, 50);
    anchorR = interp1(x_right_tail, y_AB_right, 450); % Anchored to the inner edge of the view
    refR    = anchorR * exp(-p_plus * (x_refR - 450));
    semilogy(x_refR, refR, 'k--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('slope -p_+ = -%.1f', p_plus));
        
    title('Extreme Right Tail (Log-Scale)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Log-Density', 'FontSize', 12);
    legend('Location', 'northeast', 'FontSize', 10);
    set(gca, 'YScale', 'log');
    xlim([450, 500]);
    hold off;
end