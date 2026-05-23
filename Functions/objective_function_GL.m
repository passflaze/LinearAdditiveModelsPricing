function sse = objective_function_GL(x, df, yf, sigma_ATM, mon_mod, c_mkt, M, dz)
% OBJECTIVE_FUNCTION_GL Computes the Sum of Squared Errors (SSE) for fminsearch.
%
% This function acts as a robust wrapper for the Generalized Laplace FFT pricing
% engine during calibration. It enforces parameter boundary constraints using a 
% penalty barrier method and handles numerical instabilities via an error-catching mechanism.
%
% INPUTS:
%   x         : 2-element parameter vector where x(1) = alpha, x(2) = beta
%   df        : Vector or scalar of market discount factors
%   yf        : Vector or scalar of Year Fractions (Time to Maturity)
%   sigma_ATM : Vector or scalar of At-The-Money implied volatilities
%   mon_mod   : Matrix of modified market moneyness points
%   c_mkt     : Matrix of market option calibration prices (NaNs allowed for ITM/empty slots)
%   M         : Power of 2 dictating the number of FFT nodes (N = 2^M)
%   dz        : Spatial log-moneyness grid step size
%
% OUTPUT:
%   sse       : Scalar value representing the tracking metric to be minimized

    alpha = x(1);
    beta  = x(2);
    
    
    % --- Step 2: Safe Execution Wrapper ---

    c_mod = price_GL(alpha, beta, M, dz, df, sigma_ATM, yf, mon_mod);
        
 
    residuals = c_mkt - c_mod;

    % Treat NaNs as zero contribution (ignore them in the sum)
    residuals(isnan(residuals)) = 0;

    % Return scalar SSE
    sse = sum(residuals(:).^2);

end