function exact_cdf = cdf_MA(pt_plus, pt_minus, ps_plus, ps_minus, drift, x_grid)
%CDF_MA Computes and plots the analytical CDF of a Minimal Additive increment.
%   EXACT_CDF = CDF_MA(PT_PLUS, PT_MINUS, PS_PLUS, PS_MINUS, DRIFT, X_GRID)
%   calculates the exact closed-form Cumulative Distribution Function based on 
%   the analytical proposition for the MA model. It handles the continuous 
%   exponential tails and the discrete point mass at the drift location.
%
%   Inputs:
%       pt_plus  - Right-tail decay parameter at terminal time t (pt_plus > 0)
%       pt_minus - Left-tail decay parameter at terminal time t (pt_minus > 0)
%       ps_plus  - Right-tail decay parameter at initial time s (ps_plus > 0)
%       ps_minus - Left-tail decay parameter at initial time s (ps_minus > 0)
%       drift    - The deterministic martingale correction (Delta_mu)
%       x_grid   - Array of spatial values where the CDF is evaluated
%
%   Outputs:
%       exact_cdf - Array containing the evaluated CDF probabilities

    % =========================================================================
    % STEP 1: COMPUTE ANALYTICAL CONSTANTS
    % =========================================================================
    % Point mass size 'c' (probability of zero jumps)
    c = (pt_plus * pt_minus) / (ps_plus * ps_minus);
    
    % Scaling coefficients for the continuous density branches
    denominator = pt_plus + pt_minus;
    A_plus = c * ((ps_plus - pt_plus) * (ps_minus + pt_plus)) / denominator;
    A_minus  = c * ((ps_minus - pt_minus) * (ps_plus + pt_minus)) / denominator;

    % =========================================================================
    % STEP 2: PIECEWISE CDF EVALUATION
    % =========================================================================
    % Initialize the output array
    exact_cdf = zeros(size(x_grid));
    
    % Create boolean masks to separate the domain at the drift point
    left_mask = x_grid < drift;
    right_mask = x_grid >= drift;
    
    % Evaluate the left branch (x < Delta_mu)
    exact_cdf(left_mask) = (A_minus / pt_minus) * exp(pt_minus * (x_grid(left_mask) - drift));
    
    % Evaluate the right branch (x >= Delta_mu)
    % Note the addition of 'c' representing the instantaneous vertical jump
    exact_cdf(right_mask) = (A_minus / pt_minus) + c + ...
        (A_plus / pt_plus) * (1 - exp(-pt_plus * (x_grid(right_mask) - drift)));

    % =========================================================================
    % STEP 3: VISUALIZATION
    % =========================================================================
    figure('Name', 'Analytical MA CDF');
    
    % Plot the strictly continuous segments, avoiding a slanted line at the jump
    plot(x_grid(left_mask), exact_cdf(left_mask), 'b-', 'LineWidth', 2);
    hold on;
    plot(x_grid(right_mask), exact_cdf(right_mask), 'b-', 'LineWidth', 2);
    
    % Highlight the discontinuous point mass (jump) at the drift
    left_limit = A_minus / pt_minus;
    right_value = left_limit + c;
    plot([drift, drift], [left_limit, right_value], 'r--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Point Mass Jump (c = %.4f)', c));
    
    % Formatting
    title('Closed-Form CDF of the MA Increment');
    xlabel('Spatial Grid (x)');
    ylabel('Cumulative Probability (F_{s,t}(x))');
    ylim([-0.05, 1.05]);
    grid on;
    legend('Location', 'northwest');
    hold off;

end