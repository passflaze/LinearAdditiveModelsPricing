function exact_cdf = cdf_MA(pt_plus, pt_minus, ps_plus, ps_minus, drift, x_grid)
% CDF_MA  Closed-form CDF of a Minimal Additive increment.
%   Handles exponential tails and the discrete point mass at drift (Delta_mu).
%   Also draws a figure with the jump discontinuity highlighted.
%
% INPUTS:
%   pt_plus  - right-tail decay at time t  (> 0)
%   pt_minus - left-tail  decay at time t  (> 0)
%   ps_plus  - right-tail decay at time s  (> 0)
%   ps_minus - left-tail  decay at time s  (> 0)
%   drift    - deterministic drift Delta_mu
%   x_grid   - evaluation points (vector)
% OUTPUT:
%   exact_cdf - CDF probabilities, same size as x_grid

    % Point mass probability c = ratio of normalizing constants
    c = (pt_plus * pt_minus) / (ps_plus * ps_minus);

    % Continuous branch scaling coefficients
    denominator = pt_plus + pt_minus;
    A_plus  = c * ((ps_plus  - pt_plus)  * (ps_minus + pt_plus))  / denominator;
    A_minus = c * ((ps_minus - pt_minus) * (ps_plus  + pt_minus)) / denominator;

    exact_cdf = zeros(size(x_grid));

    left_mask  = x_grid < drift;
    right_mask = x_grid >= drift;

    exact_cdf(left_mask) = (A_minus / pt_minus) * exp(pt_minus * (x_grid(left_mask) - drift));

    % Right branch: continuous part plus the point mass atom at drift
    exact_cdf(right_mask) = (A_minus / pt_minus) + c + ...
        (A_plus / pt_plus) * (1 - exp(-pt_plus * (x_grid(right_mask) - drift)));

    figure('Name', 'Analytical MA CDF');

    % Plot the two continuous segments separately to avoid a spurious slanted
    % line across the jump discontinuity at drift.
    plot(x_grid(left_mask),  exact_cdf(left_mask),  'b-', 'LineWidth', 2);
    hold on;
    plot(x_grid(right_mask), exact_cdf(right_mask), 'b-', 'LineWidth', 2);

    left_limit  = A_minus / pt_minus;
    right_value = left_limit + c;
    plot([drift, drift], [left_limit, right_value], 'r--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Point Mass Jump (c = %.4f)', c));

    title('Closed-Form CDF of the MA Increment');
    xlabel('Spatial Grid (x)');
    ylabel('Cumulative Probability (F_{s,t}(x))');
    ylim([-0.05, 1.05]);
    grid on;
    legend('Location', 'northwest');
    hold off;

end