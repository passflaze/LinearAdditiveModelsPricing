function [cdf_exact, x_grid] = cdf_MA_increments(pt_plus, pt_minus, ps_plus, ps_minus, delta_mu)
% CDF_MA_INCREMENTS  Closed-form CDF of the MA finite-activity increment.
%   Generates a high-resolution spatial grid and evaluates the piecewise CDF.
%   Also draws a figure highlighting the point-mass jump at delta_mu.
%
% INPUTS:
%   pt_plus  - right-tail decay at t  (> 0)
%   pt_minus - left-tail  decay at t  (> 0)
%   ps_plus  - right-tail decay at s  (> 0)
%   ps_minus - left-tail  decay at s  (> 0)
%   delta_mu - deterministic drift of the increment
% OUTPUTS:
%   cdf_exact - (column vector) CDF values on x_grid
%   x_grid    - (column vector) spatial grid, 10^7 points on [-1000, 1000]

    % High-resolution spatial grid
    x_grid = linspace(-1000, 1000, 10000000)';

    % Point mass probability c = ratio of normalizing constants
    c = (pt_plus * pt_minus) / (ps_plus * ps_minus);

    % Continuous branch scaling coefficients
    denominator = pt_plus + pt_minus;
    A_plus  = c * ((ps_plus  - pt_plus)  * (ps_minus + pt_plus))  / denominator;
    A_minus = c * ((ps_minus - pt_minus) * (ps_plus  + pt_minus)) / denominator;

    cdf_exact = zeros(size(x_grid));

    idx_left = x_grid < delta_mu;
    cdf_exact(idx_left) = (A_minus / pt_minus) .* exp(pt_minus .* (x_grid(idx_left) - delta_mu));

    idx_right = x_grid >= delta_mu;
    term_right = (A_plus / pt_plus) .* (1 - exp(-pt_plus .* (x_grid(idx_right) - delta_mu)));
    cdf_exact(idx_right) = (A_minus / pt_minus) + c + term_right;

    figure('Name', 'Analytical CDF of MA Increment', 'Color', 'w');
    plot(x_grid, cdf_exact, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical CDF');
    hold on;

    xline(delta_mu, 'k:', 'LineWidth', 1.5, 'DisplayName', '\Delta\mu Location');

    % Vertical segment marking the point-mass jump at delta_mu
    y_bottom = A_minus / pt_minus;
    y_top    = y_bottom + c;
    plot([delta_mu, delta_mu], [y_bottom, y_top], 'r-', 'LineWidth', 3, 'DisplayName', 'Point Mass (c)');

    title('Analytical CDF of the Finite-Activity MA Increment');
    xlabel('Spatial Grid (x)');
    ylabel('Cumulative Probability F(x)');
    xlim([-100 100]);
    legend('Location', 'best');
    grid on;
    hold off;

end