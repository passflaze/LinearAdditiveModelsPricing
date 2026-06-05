function X_samples = simulate_from_cdf(cdf_clean, x_grid_fine, spline, N_sim)
%SIMULATE_FROM_CDF Inverse-CDF sampling from a pre-computed numerical CDF.
%   X_SAMPLES = SIMULATE_FROM_CDF(CDF_CLEAN, X_GRID_FINE, SPLINE, N_SIM) draws
%   N_sim independent samples by numerically inverting a CDF that has already
%   been reconstructed (e.g. via lewis_FFT_digital / tail_adjustment). Two
%   inversion schemes are available:
%     * spline = true  : PCHIP inversion. The (CDF, x) pairs are made strictly
%                         increasing with unique(), then U ~ U(0,1) is mapped
%                         through interp1(cdf, x, U, 'pchip', 'extrap'). This is
%                         the path used throughout the project.
%     * spline = false : piecewise-linear inversion (Glasserman & Liu, 2010):
%                         bracket each U between adjacent CDF nodes via
%                         discretize() and invert the local linear segment.
%
%   Inputs:
%       cdf_clean   - CDF values (monotone, 0 -> 1) on x_grid_fine
%       x_grid_fine - spatial grid corresponding to cdf_clean
%       spline      - (logical) true for PCHIP inversion, false for linear
%       N_sim       - number of independent samples to generate
%
%   Outputs:
%       X_samples - column vector of N_sim simulated increments
%
%   References:
%       Glasserman and Liu (2010) linear interpolation scheme (spline = false).


    if spline

        [cdf_clean,unique_idx] = unique(cdf_clean,'stable');
        x_grid_fine = x_grid_fine(unique_idx);
        U = rand(N_sim, 1);
        X_samples = interp1(cdf_clean,x_grid_fine,U,'pchip','extrap');

    else
        % Force exact bounds so every U is bracketed.
        cdf_clean(1) = 0;
        cdf_clean(end) = 1;

        U = rand(N_sim, 1);

        % Bracket each U between adjacent CDF nodes; clamp U = 1 to the last bin.
        j = discretize(U, cdf_clean) + 1;
        j(j > length(cdf_clean)) = length(cdf_clean);

        % Bounding CDF values and x-coordinates of each bracket.
        P_j = cdf_clean(j);
        P_j_minus_1 = cdf_clean(j - 1);

        X_j = x_grid_fine(j);
        X_j_minus_1 = x_grid_fine(j - 1);

        gamma = X_j - X_j_minus_1;

        % delta_P > 0 because the CDF is strictly monotonic after tail_adjustment.
        delta_P = P_j - P_j_minus_1;

        c_1 = gamma ./ delta_P;
        c_0 = (X_j .* P_j_minus_1 - X_j_minus_1 .* P_j) ./ delta_P;

        % Local linear inversion: X = c_0 + c_1 * U.
        X_samples = c_0 + c_1 .* U;
    end

end