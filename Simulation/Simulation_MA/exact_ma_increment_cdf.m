function [cdf_exact, x_grid] = exact_ma_increment_cdf(pt_plus, pt_minus, ps_plus, ps_minus, delta_mu)
%EXACT_MA_INCREMENT_CDF Computes the analytical CDF of the MA increment.
%   [CDF_EXACT, X_GRID] = EXACT_MA_INCREMENT_CDF(PT_PLUS, PT_MINUS, ...
%   PS_PLUS, PS_MINUS, DELTA_MU) generates a spatial grid from -40 to 40
%   and calculates the piecewise exact Cumulative Distribution Function 
%   based on the finite-activity jump distribution.

    % Generazione della griglia spaziale (100.000 punti per alta risoluzione)
    x_grid = linspace(-1000, 1000, 10000000)';

    % Point mass size 'c' (probability of zero jumps)
    c = (pt_plus * pt_minus) / (ps_plus * ps_minus);
    
    % Scaling coefficients for the continuous density branches
    denominator = pt_plus + pt_minus;
    A_plus = c * ((ps_plus - pt_plus) * (ps_minus + pt_plus)) / denominator;
    A_minus  = c * ((ps_minus - pt_minus) * (ps_plus + pt_minus)) / denominator;

    % Inizializzazione del vettore CDF
    cdf_exact = zeros(size(x_grid));

    % Vettorizzazione logica per x < \Delta\mu
    idx_left = x_grid < delta_mu;
    cdf_exact(idx_left) = (A_minus / pt_minus) .* exp(pt_minus .* (x_grid(idx_left) - delta_mu));

    % Vettorizzazione logica per x >= \Delta\mu
    idx_right = x_grid >= delta_mu;
    term_right = (A_plus / pt_plus) .* (1 - exp(-pt_plus .* (x_grid(idx_right) - delta_mu)));
    cdf_exact(idx_right) = (A_minus / pt_minus) + c + term_right;

    figure('Name', 'Analytical CDF of MA Increment', 'Color', 'w');
    plot(x_grid, cdf_exact, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical CDF');
    hold on;
    
    % Evidenziamo il punto di salto (\Delta\mu)
    xline(delta_mu, 'k:', 'LineWidth', 1.5, 'DisplayName', '\Delta\mu Location');
    
    % Evidenziamo il gradino causato dalla massa discreta 'c'
    y_bottom = A_minus / pt_minus;
    y_top = y_bottom + c;
    plot([delta_mu, delta_mu], [y_bottom, y_top], 'r-', 'LineWidth', 3, 'DisplayName', 'Point Mass (c)');
    
    title('Analytical CDF of the Finite-Activity MA Increment');
    xlabel('Spatial Grid (x)');
    ylabel('Cumulative Probability F(x)');
    xlim([-100 100]); % Limita la visuale per vedere meglio il salto (modificabile)
    legend('Location', 'best');
    grid on;
    hold off;

end