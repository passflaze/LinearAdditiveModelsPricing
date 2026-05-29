% Initialization
clc;
clear;
close all;

%% 1. Global Parameters and Evaluation Grid
% Setting universal tail decay parameters for consistency across models
p_plus  = 1.5;
p_minus = 0.9;

% Unified sampling grid
x_grid = linspace(-30, 30, 500); 

%% 2. Asymmetric Laplace Distribution (MA Model)
% Parameter mapping
alpha_MA = p_minus;
beta_MA  = p_plus;

% Analytical components
gamma_MA = @(a, b) (1/a - 1/b);
C_MA     = @(a, b) (1/b + 1/a)^(-1);
A_MA     = @(a, b) (a - b) / 2;
B_MA     = @(a, b) (a + b) / 2;

% PDF definition (fully vectorized)
pdf_MA = @(a, b, x) C_MA(a,b) .* exp(A_MA(a,b).*(x - gamma_MA(a,b)) - B_MA(a,b).*abs(x - gamma_MA(a,b)));

% Evaluation
y_MA = pdf_MA(alpha_MA, beta_MA, x_grid);

%% AB PDF via NIG Model Simulation
% Centralized parameters (assumed from workspace setup)
p_plus  = 1.5;
p_minus = 0.9;

kappa_AB = 1 / (p_plus * p_minus);
eta_AB   = ((1/kappa_AB)-(p_minus^2)) / (2 * p_minus);

% NIG structural parameter mapping
mu      = eta_AB;
gamma_AB   = -eta_AB;
mu_IG1  = 1;
lam_IG1 = 1 / kappa_AB;

alpha_nig = sqrt(gamma_AB^2 + lam_IG1);
beta_nig = gamma_AB;
delta_nig = sqrt(lam_IG1);
mu_nig = mu;

y_AB = nigpdf(x_grid, alpha_nig, beta_nig, mu_nig, delta_nig);

% Grafico della PDF continua
% 1. Correctly open the figure window first
figure('Color', 'w', 'Position', [100, 100, 900, 600]);

% 2. Activate the hold state on the current axis
hold on;
grid on;
box on;

% 3. Plot the empirical KDE curve
plot(x_grid, y_AB, 'LineWidth', 2.5, 'Color', [0.4940 0.1840 0.5560], ...
     'DisplayName', 'Empirical PDF (KDE)');

% % 4. Plot the analytical Gil-Pelaez curve (Fixed invalid RGB values to a professional dark orange/red)
% plot(x_grid, y_AB, '--','LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980], ...
%      'DisplayName', 'Analytical PDF (Gil-Pelaez)');

% 5. Add formatting, typography, and layout
title('AB Distribution: Empirical KDE vs. Analytical PDF', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Density', 'FontSize', 12);

% 6. Add a professional legend using LaTeX formatting
legend('Location', 'northeast', 'FontSize', 11, 'Interpreter', 'latex');

% 7. Enforce exact range boundaries and release hold
xlim([min(x_grid), max(x_grid)]);
hold off;


%% pdfnig vs simulazione

% Simulation size
numSim = 1e6; 

% Sample the Subordinator (Time-Change) from the Inverse Gaussian
W = random('InverseGaussian', mu_IG1, lam_IG1, [numSim, 1]);

% Sample the standard normal innovations (randn is much faster than random)
Z = randn(numSim, 1);

% Construct the AB/NIG process (Variance-Mean Mixture)
y_sim_AB = mu + W .* gamma_AB + sqrt(W) .* Z;

% Compute empirical PDF via Kernel Density Estimation
[f_emp_AB, ~] = ksdensity(y_sim_AB, x_grid);

% Plotting (Linear Scale)
figure('Color', 'w', 'Position', [150, 150, 900, 600]);
hold on; 
grid on; 
box on;

% Plot Empirical PDF (Simulated)
plot(x_grid, f_emp_AB, 'LineWidth', 4, 'Color', [0.4940 0.1840 0.5560], ...
     'DisplayName', 'Empirical PDF (Monte Carlo 1M)');

% Plot Analytical PDF (From custom nigpdf function)
plot(x_grid, y_AB, '--', 'LineWidth', 2.5, 'Color', [0.9290 0.6940 0.1250], ...
     'DisplayName', 'Analytical PDF (nigpdf)');

title('Additive Bachelier: Monte Carlo Simulation vs Analytical PDF', 'FontSize', 15, 'FontWeight', 'bold');
xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Density', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');

xlim([-8, 8]); 
ylim([0, max(y_AB) * 1.1]);
hold off;

%% 4. Generalized Logistic Distribution (GL)
% Parameter mapping
alpha_GL = p_minus;
beta_GL  = p_plus;

% Translation parameter (using MATLAB's built-in digamma function 'psi')
gamma_GL = psi(beta_GL) - psi(alpha_GL); 

% Normalization constant
C_GL = gamma(alpha_GL + beta_GL) / (gamma(alpha_GL) * gamma(beta_GL));

% PDF definition (fully vectorized)
pdf_GL = @(a, b, g, x) C_GL .* (exp(a .* (x - g)) ./ (1 + exp(x - g)).^(a + b));

% Evaluation
y_GL = pdf_GL(alpha_GL, beta_GL, gamma_GL, x_grid);


%% 4.5. Final Comparison Plot linear

% 1. Define Normal Distribution parameters (Standard Normal as benchmark)
mu_norm = 0;
sigma_norm = 1;

% 2. Calculate Normal PDF (Explicit formula to avoid Statistics Toolbox dependency)
y_norm = (1 / (sigma_norm * sqrt(2*pi))) .* exp(-0.5 .* ((x_grid - mu_norm) ./ sigma_norm).^2);

% 3. Create the figure
figure('Color', 'w', 'Position', [100, 100, 900, 600]);
hold on; 
grid on; 
box on;

% Plotting distributions with standard MATLAB professional color palette
plot(x_grid, y_MA, 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Asymmetric Laplace (MA)');
plot(x_grid, y_AB, 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB Distribution');
plot(x_grid, y_GL, 'LineWidth', 2.5, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'Generalized Logistic (GL)');

% Adding the Normal curve (Dashed black line for comparison)
plot(x_grid, y_norm, '--k', 'LineWidth', 2, 'DisplayName', 'Normal Distribution Benchmark');

% Formatting and aesthetics
title('Probability Density Functions Comparison (Linear Scale)', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Density', 'FontSize', 14);

% Legend layout
legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');

% Axis limits (Updated to ensure the Normal peak is included if it's the highest)
xlim([min(x_grid), max(x_grid)]);
ylim([0, max([y_MA, y_AB, y_GL, y_norm]) * 1.1]);

hold off;




%% 5. Final Comparison Plot
figure('Color', 'w', 'Position', [100, 100, 900, 600]);
hold on; 
grid on; 
box on;

% Plotting distributions with standard MATLAB professional color palette
semilogy(x_grid, y_MA, 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Asymmetric Laplace (MA)');
semilogy(x_grid, y_AB, 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB Distribution');
semilogy(x_grid, y_GL, 'LineWidth', 2.5, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'Generalized Logistic (GL)');
set(gca, 'YScale', 'log');

% Formatting and aesthetics
title('Probability Density Functions Comparison semilog-scale', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Density', 'FontSize', 14);

% Legend layout
legend('Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');

% Axis limits
xlim([min(x_grid), max(x_grid)]);
ylim([1e-10, max([y_MA, y_AB, y_GL]) * 1.1]);

hold off;
%% 6. Tail Analysis (Zoom with Semi-Logarithmic Scale)
% In quantitative finance, tails are best visualized in log-scale 
% to verify the linear asymptotic exponential decay.

figure('Color', 'w', 'Position', [150, 150, 1000, 450]);

% --- Subplot 1: Left Tail (Zoom) ---
subplot(1, 2, 1);
hold on; grid on; box on;
semilogy(x_grid, y_MA, 'LineWidth', 2, 'Color', [0 0.4470 0.7410], 'DisplayName', 'MA');
semilogy(x_grid, y_AB, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB');
semilogy(x_grid, y_GL, 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'GL');

title('Left Tail (Log-Scale)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Log-Density', 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 10);
set(gca, 'YScale', 'log');
% Zooming in on the extreme negative values
xlim([-30, -10]); 
ylim([1e-10, 1e-3]); % Setting a fixed floor for visibility

% --- Subplot 2: Right Tail (Zoom) ---
subplot(1, 2, 2);
hold on; grid on; box on;
semilogy(x_grid, y_MA, 'LineWidth', 2, 'Color', [0 0.4470 0.7410], 'DisplayName', 'MA');
semilogy(x_grid, y_AB, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'AB');
semilogy(x_grid, y_GL, 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'GL');

title('Right Tail (Log-Scale)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Log-Density', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'YScale', 'log');

% Zooming in on the extreme positive values
xlim([10, 30]); 
ylim([1e-10, 1e-3]);