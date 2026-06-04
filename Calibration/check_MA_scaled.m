% CHECK_MA_SCALED  Standalone demo: MA density at calibrated (alpha,beta) vs the
%   same rescaled by 1/10, confirming the gauge invariance of the scaled density.

set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');
set(groot, 'defaultAxesXColor', 'k');
set(groot, 'defaultAxesYColor', 'k');
set(groot, 'defaultTextColor', 'k');

x_grid = linspace(-100, 100, 1000);

% Calibrated MA parameters and ATM vol (hardcoded for the demo).
alphaMA1 = 1.000000000000000;
betaMA1  = 0.982156170515665;
sigmaATM = 14.740269016504566;

% Distribution 1: scaled density f(sigma*x) = (1/sigma) f(x/sigma).
sigmaMA1 = sigmaATM / I0_MA([alphaMA1; betaMA1]);
y_MA1 = (1 / sigmaMA1) * pdf_MA([alphaMA1; betaMA1], 1, x_grid / sigmaMA1);

% Distribution 2: parameters scaled down by 10.
alphaMA2 = alphaMA1 / 10;
betaMA2  = betaMA1 / 10;
sigmaMA2 = sigmaATM / I0_MA([alphaMA2; betaMA2]);
y_MA2 = (1 / sigmaMA2) * pdf_MA([alphaMA2; betaMA2], 1, x_grid / sigmaMA2);

figure('Name', 'MA Distributions Scale Comparison', 'Position', [150, 150, 900, 500]);
hold on; grid on; box on;

plot(x_grid, y_MA1, '-', 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', sprintf('MA (\\alpha=%.1f, \\beta\\approx%.2f)', alphaMA1, betaMA1));
plot(x_grid, y_MA2, '--', 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', sprintf('MA (\\alpha/10, \\beta/10)'));

title('Comparison of Scaled Asymmetric Laplace (MA) Distributions', 'FontSize', 15, 'FontWeight', 'bold');
xlabel('Modified Moneyness ($\zeta$)', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Scaled Density', 'FontSize', 13);
legend('Location', 'best', 'FontSize', 11);
xlim([min(x_grid), max(x_grid)]);

hold off;