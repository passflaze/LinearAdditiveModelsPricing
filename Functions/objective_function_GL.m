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
    
    % --- Step 1: Execute Pricing Engine ---
    c_mod = price_GL(alpha, beta, M, dz, df, sigma_ATM, yf, mon_mod);
    
    % --- Step 2: STRICT PENALTY CHECK ---
    % Se c_mod è NaN (fallimento del modello) MA c_mkt è un numero valido, scatta la sanzione.
    if any(isnan(c_mod(:)) & ~isnan(c_mkt(:)))
        fprintf('[PENALTY] Il modello ha restituito NaN su un target valido di mercato. Alpha: %.4f, Beta: %.4f\n', alpha, beta);
        sse = 1e10 + 1e5 * ((alpha - 1)^2 + (beta - 1)^2);
        return;
    end
    
    % --- Step 3: SAFE SSE CALCULATION ---
    % Sottraiamo le matrici. Se in una cella c'è un NaN (es. l'opzione di mercato non c'era),
    % la differenza rimarrà semplicemente NaN.
    diff = c_mkt - c_mod;
    
    % Creiamo una maschera che seleziona SOLO i punti in cui la differenza è un numero vero.
    valid_idx = ~isnan(diff);
    
    % Controllo anti-collasso (evita il finto "SSE = 0" che inganna fmincon)
    if sum(valid_idx(:)) == 0 || all(c_mod(valid_idx) == 0)
        fprintf('[PENALTY] La superficie del modello è collassata a zero. Alpha: %.4f, Beta: %.4f\n', alpha, beta);
        sse = 1e10 + 1e5 * ((alpha - 1)^2 + (beta - 1)^2);
        return;
    end
    
    % Estraiamo solo i residui validi (diventa un vettore) ed eleviamo al quadrato
    residuals = diff(valid_idx);
    sse = sum(residuals.^2);
    
end