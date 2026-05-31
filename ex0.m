% Initialization
addpath('Distributions');
clc;
clear;
close all;

%% 1. Global Parameters and Evaluation Grid
p_plus  = 1.5;
p_minus = 0.9;

x_grid = linspace(-30, 30, 500);

%% 2. MA PDF (Asymmetric Laplace)
% sigma*sqrt(t) = 1 for unit-time evaluation
y_MA = pdf_MA([p_minus; p_plus], 1, x_grid);

%% 3. AB PDF (NIG representation)
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

%% 4. GL PDF (Generalized Logistic)
y_GL = pdf_GL([p_minus; p_plus], x_grid);

%% AB Monte Carlo Simulation vs Analytical PDF
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

%% PDF Comparison (Linear Scale)
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
xlim([min(x_grid), max(x_grid)]);
ylim([0, max([y_MA, y_AB, y_GL, y_norm]) * 1.1]);
hold off;

%% PDF Comparison (Semi-Log Scale)
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
xlim([min(x_grid), max(x_grid)]);
ylim([1e-10, max([y_MA, y_AB, y_GL]) * 1.1]);
hold off;

%% Tail Analysis (Zoom, Semi-Log Scale)
figure('Color', 'w', 'Position', [150, 150, 1000, 450]);

subplot(1, 2, 1);
hold on; grid on; box on;
semilogy(x_grid, y_MA, 'LineWidth', 2, 'Color', [0 0.4470 0.7410],     'DisplayName', 'MA');
semilogy(x_grid, y_AB, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB');
semilogy(x_grid, y_GL, 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'GL');
% Reference asymptotic slope (point 0a): left tail ~ exp(+p_minus * x), a straight
% line of slope p_minus on a semilog plot. Anchored to the AB pdf at x = -10; the
% three model tails should run parallel to it.
x_refL  = linspace(-30, -10, 50);
anchorL = interp1(x_grid, y_AB, -10);
refL    = anchorL * exp(p_minus * (x_refL + 10));
semilogy(x_refL, refL, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('slope p_- = %.1f', p_minus));
title('Left Tail (Log-Scale)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Log-Density', 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 10);
set(gca, 'YScale', 'log');
xlim([-30, -10]);
ylim([1e-10, 1e-3]);

subplot(1, 2, 2);
hold on; grid on; box on;
semilogy(x_grid, y_MA, 'LineWidth', 2, 'Color', [0 0.4470 0.7410],     'DisplayName', 'MA');
semilogy(x_grid, y_AB, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB');
semilogy(x_grid, y_GL, 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'GL');
% Reference asymptotic slope (point 0a): right tail ~ exp(-p_plus * x), a straight
% line of slope -p_plus on a semilog plot. Anchored to the AB pdf at x = +10; the
% three model tails should run parallel to it.
x_refR  = linspace(10, 30, 50);
anchorR = interp1(x_grid, y_AB, 10);
refR    = anchorR * exp(-p_plus * (x_refR - 10));
semilogy(x_refR, refR, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('slope -p_+ = -%.1f', p_plus));
title('Right Tail (Log-Scale)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Log-Density', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'YScale', 'log');
xlim([10, 30]);
ylim([1e-10, 1e-3]);
