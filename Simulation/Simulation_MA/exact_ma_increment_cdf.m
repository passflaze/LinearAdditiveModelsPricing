function [cdf_exact, x_grid] = exact_ma_increment_cdf(pt_plus, pt_minus, ps_plus, ps_minus, delta_mu, doPlot)
%EXACT_MA_INCREMENT_CDF Computes the analytical CDF of the MA increment.
%   [CDF_EXACT, X_GRID] = EXACT_MA_INCREMENT_CDF(PT_PLUS, PT_MINUS, PS_PLUS,
%   PS_MINUS, DELTA_MU) evaluates, on a fine spatial grid, the piecewise exact
%   Cumulative Distribution Function of the finite-activity MA increment over
%   [s, t]. The distribution has a point mass 'c' at DELTA_MU (the no-jump
%   atom) plus two exponentially decaying continuous branches on either side.
%
%   INPUTS:
%       pt_plus, pt_minus : right/left tail decay parameters at terminal time t
%       ps_plus, ps_minus : right/left tail decay parameters at initial time s
%       delta_mu          : deterministic drift of the increment (atom location)
%       doPlot            : (optional, default false) plot the CDF, the atom
%                           location and the discrete-mass step
%
%   OUTPUTS:
%       cdf_exact         : exact CDF values on x_grid
%       x_grid            : spatial grid (column vector)

    if nargin < 6 || isempty(doPlot)
        doPlot = false;
    end

    % Spatial grid: +-200 already pushes 1 - CDF below ~1e-30 for the relevant
    % tail decays, so 2e5 points give ample resolution without the memory and
    % runtime cost of an oversized grid (the previous 1e7-point grid made the
    % downstream interp1/integral verification needlessly slow).
    x_grid = linspace(-200, 200, 2e5)';

    % Point mass size 'c' (probability of zero jumps)
    c = (pt_plus * pt_minus) / (ps_plus * ps_minus);

    % Scaling coefficients for the continuous density branches
    denominator = pt_plus + pt_minus;
    A_plus = c * ((ps_plus - pt_plus) * (ps_minus + pt_plus)) / denominator;
    A_minus  = c * ((ps_minus - pt_minus) * (ps_plus + pt_minus)) / denominator;

    % CDF vector initialization
    cdf_exact = zeros(size(x_grid));

    % Left branch (x < delta_mu): only left jumps contribute
    idx_left = x_grid < delta_mu;
    cdf_exact(idx_left) = (A_minus / pt_minus) .* exp(pt_minus .* (x_grid(idx_left) - delta_mu));

    % Right branch (x >= delta_mu): left-tail mass + atom + right-jump mass
    idx_right = x_grid >= delta_mu;
    term_right = (A_plus / pt_plus) .* (1 - exp(-pt_plus .* (x_grid(idx_right) - delta_mu)));
    cdf_exact(idx_right) = (A_minus / pt_minus) + c + term_right;

    if doPlot
        figure('Name', 'Analytical CDF of MA Increment', 'Color', 'w');
        plot(x_grid, cdf_exact, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical CDF');
        hold on;

        % Mark the atom location (delta_mu)
        xline(delta_mu, 'k:', 'LineWidth', 1.5, 'DisplayName', '\Delta\mu Location');

        % Mark the jump caused by the discrete mass 'c'
        y_bottom = A_minus / pt_minus;
        y_top = y_bottom + c;
        plot([delta_mu, delta_mu], [y_bottom, y_top], 'r-', 'LineWidth', 3, 'DisplayName', 'Point Mass (c)');

        title('Analytical CDF of the Finite-Activity MA Increment');
        xlabel('Spatial Grid (x)');
        ylabel('Cumulative Probability F(x)');
        xlim([-100 100]); % Zoom in to see the jump (adjustable)
        legend('Location', 'best');
        grid on;
        hold off;
    end

end