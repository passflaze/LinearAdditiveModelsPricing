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
    % =========================================================================
    % STEP 1: GENERATE THE CLEAN NUMERICAL CDF
    % =========================================================================
    % Call the robust FFT inversion routine (which includes tail adjustment)
        
        
        cdf_clean(1) = 0;
        cdf_clean(end) = 1;
    
        % =========================================================================
        % STEP 2: GENERATE UNIFORM RANDOM VARIABLES
        % =========================================================================
        % Generate N_sim standard uniform random variables U ~ U(0,1)
        U = rand(N_sim, 1);
    
        % =========================================================================
        % STEP 3: FIND NEAREST NEIGHBORHOOD (INDEX MAPPING)
        % =========================================================================
        j = discretize(U, cdf_clean) + 1;
        
        % Handle edge case: if U = 1 exactly, it might map outside, clamp it to end
        j(j > length(cdf_clean)) = length(cdf_clean);
    
        % =========================================================================
        % STEP 4: EXTRACT LOCAL GRID POINTS
        % =========================================================================
        % For each simulated U, extract the bounding CDF values and X coordinates
        P_j = cdf_clean(j);
        P_j_minus_1 = cdf_clean(j - 1);
        
        X_j = x_grid_fine(j);
        X_j_minus_1 = x_grid_fine(j - 1);
        
        % Calculate the spatial step size (gamma in the paper)
        gamma = X_j - X_j_minus_1;
    
        % =========================================================================
        % STEP 5: COMPUTE INTERPOLATION COEFFICIENTS
        % =========================================================================
        % Note: Since our CDF is strictly monotonic after tail_adjustment, 
        % the denominator (P_j - P_j_minus_1) is strictly > 0, avoiding DivByZero.
        
        delta_P = P_j - P_j_minus_1;
        
        c_1 = gamma ./ delta_P;
        c_0 = (X_j .* P_j_minus_1 - X_j_minus_1 .* P_j) ./ delta_P;
    
        % =========================================================================
        % STEP 6: COMPUTE FINAL SAMPLES
        % =========================================================================
        % Linear inversion formula: X = c_0 + c_1 * U
        X_samples = c_0 + c_1 .* U;
    end

end