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
    
    % --- Step 1: Boundary Constraint Verification ---
    if alpha <= 1e-3 || beta <= 1e-3 || alpha > 50 || beta > 50
        fprintf('[PENALTY] Constraints violated. Alpha: %.4f, Beta: %.4f (Must be between 0.001 and 50).\n', alpha, beta);
        sse = 1e10; 
        return; 
    end
    
    % --- Step 2: Safe Execution Wrapper ---
    try
        c_mod = price_GL(alpha, beta, M, dz, df, sigma_ATM, yf, mon_mod);
        
        % Check if pricing engine returned completely invalid matrices
        if any(isnan(c_mod(:))) && ~all(isnan(mon_mod(:)))
            fprintf('[PENALTY] Model surface contains NaNs at Alpha: %.4f, Beta: %.4f.\n', alpha, beta);
            sse = 1e10;
            return;
        end
        
    catch ME
        % Inform exactly which component crashed and print the error stack message
        fprintf('[CRASH] Execution failed inside pricing engine at Alpha: %.4f, Beta: %.4f.\n', alpha, beta);
        fprintf('        Error Identifier: %s\n', ME.identifier);
        fprintf('        Error Message   : %s\n', ME.message);
        sse = 1e10;
        return;
    end
    
    % --- Step 3: SSE Calculation ---
    diff = c_mkt - c_mod;
    valid_idx = ~isnan(diff);
    residuals = diff(valid_idx);
    
    sse = sum(residuals.^2);
    fprintf('[SUCCESS] Alpha: %.4f, Beta: %.4f -> Current SSE: %.6e\n', alpha, beta, sse);
end