function X_samples = simulate_from_cdf(cdf_clean, x_grid_fine, spline,N_sim)
%SIMULATE_FROM_CDF Generates random samples from the MA process via FFT inversion.
%   X_SAMPLES = SIMULATE_FROM_CDF(N_SIM, M, DX, SHIFT, PT_PLUS, PT_MINUS, 
%   PS_PLUS, PS_MINUS) computes the numerical CDF of a Minimal Additive 
%   increment and generates N_sim independent random samples by numerically 
%   inverting the CDF using vectorized linear interpolation.
%
%   Inputs:
%       N_sim     - Number of Monte Carlo paths/samples to generate
%       M         - Controls grid size for FFT, where N = 2^M
%       dx        - Step size for the spatial grid (x-domain)
%       shift     - Damping factor (alpha)
%       pt_plus   - Right-tail decay parameter at terminal time t
%       pt_minus  - Left-tail decay parameter at terminal time t
%       ps_plus   - Right-tail decay parameter at initial time s
%       ps_minus  - Left-tail decay parameter at initial time s
%
%   Outputs:
%       X_samples - A vector of N_sim simulated random increments
%
%   References:
%       Glasserman and Liu (2010) linear interpolation scheme.
    

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
        % Use 'discretize' to perform a fast O(N_sim * log(K)) binary search.
        % It returns the index j-1 such that cdf_clean(j-1) <= U < cdf_clean(j).
        % We add 1 to get the upper index 'j' to match the paper's notation.
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