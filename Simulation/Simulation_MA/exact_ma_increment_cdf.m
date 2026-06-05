function [cdf_exact, x_grid] = exact_ma_increment_cdf(pt_plus, pt_minus, ps_plus, ps_minus, delta_mu, doPlot)
%EXACT_MA_INCREMENT_CDF Computes the analytical CDF of the MA increment.
%   Builds a spatial grid and evaluates the piecewise-exact CDF of the
%   finite-activity MA increment, with a point mass at the drift delta_mu.
%
% INPUTS:
%   pt_plus  - right-tail decay parameter at terminal time t
%   pt_minus - left-tail decay parameter at terminal time t
%   ps_plus  - right-tail decay parameter at initial time s
%   ps_minus - left-tail decay parameter at initial time s
%   delta_mu - deterministic drift (location of the point mass)
%   doPlot   - (optional, default false) plot the resulting CDF
%
% OUTPUTS:
%   cdf_exact - analytical CDF values on x_grid
%   x_grid    - spatial grid the CDF is evaluated on

    % --- Optional Arguments Initialization ---
    if nargin < 6 || isempty(doPlot)
        doPlot = false;
    end

    % Spatial grid generation 
    x_grid = linspace(-1000, 1000, 10000000)';
    
    % Point mass size 'c' 
    c = (pt_plus * pt_minus) / (ps_plus * ps_minus);

    % Scaling coefficients for the continuous density branches
    denominator = pt_plus + pt_minus;
    A_plus   = c * ((ps_plus - pt_plus) * (ps_minus + pt_plus)) / denominator;
    A_minus  = c * ((ps_minus - pt_minus) * (ps_plus + pt_minus)) / denominator;
    
    % Initialization of the CDF vector
    cdf_exact = zeros(size(x_grid));
    
    % Logical vectorization for x < \Delta\mu
    idx_left = x_grid < delta_mu;
    cdf_exact(idx_left) = (A_minus / pt_minus) .* exp(pt_minus .* (x_grid(idx_left) - delta_mu));
    
    % Logical vectorization for x >= \Delta\mu
    idx_right = x_grid >= delta_mu;
    term_right = (A_plus / pt_plus) .* (1 - exp(-pt_plus .* (x_grid(idx_right) - delta_mu)));
    cdf_exact(idx_right) = (A_minus / pt_minus) + c + term_right;
    
    % --- Plotting Block ---
    if doPlot
        figure('Name', 'Analytical CDF of MA Increment', 'Color', 'w');
        plot(x_grid, cdf_exact, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical CDF');
        hold on;
        
        xline(delta_mu, 'k:', 'LineWidth', 1.5, 'DisplayName', '\Delta\mu Location');
        
        y_bottom = A_minus / pt_minus;
        y_top = y_bottom + c;
        plot([delta_mu, delta_mu], [y_bottom, y_top], 'r-', 'LineWidth', 3, 'DisplayName', 'Point Mass (c)');
        
        title('Analytical CDF of the Finite-Activity MA Increment');
        xlabel('Spatial Grid (x)');
        ylabel('Cumulative Probability F(x)');
        xlim([-100 100]); 
        legend('Location', 'best');
        grid on;
        hold off;
    end
end