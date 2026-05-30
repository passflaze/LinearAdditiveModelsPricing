function cf = cf_GL_cal(z,alpha, beta)
% CF_GL Computes the Characteristic Function of the Generalized Laplace distribution.
%
% This function evaluates the theoretical characteristic function (CF) derived from
% the ratio of Beta functions, modified with a complex exponential drift.
%
% Mathematical Formula:
%   phi(z) = [ B(alpha + i*z, beta - i*z) / B(alpha, beta) ] * exp(i * gamma_GL * z)
% 
% INPUTS:
%   alpha    : Left shape parameter (must be > 0)
%   beta     : Right shape parameter (must be > 0)
%   z        : Vector or matrix of frequencies (can be real or complex)
%   gamma_GL : Location/Drift parameter (real number)
%
% OUTPUT:
%   cf       : Vector/matrix of evaluated characteristic function values.
    
% --- Step 1: Verify Input Integrity ---
    if any(imag(alpha) ~= 0) || any(imag(beta) ~= 0)
        error('CharacteristicFunction:ComplexParameters', 'Alpha and Beta shape variables must be strictly Real inputs.');
    end

    % --- Step 2: Compute Structural Shift and Factors ---
    gamma_GL = psi(beta) - psi(alpha);
    
    gamma_factor_1 = (complex_gamma(alpha + 1i.*z) .* complex_gamma(beta - 1i.*z)) ./ complex_gamma(alpha + beta);
    gamma_factor_2 = exp(betaln(alpha, beta));
    
    % --- Step 3: Assemble Final Vectorized Profile ---
    cf = (gamma_factor_1 ./ gamma_factor_2) .* exp(1i .* gamma_GL .* z);
    nan_mask = (isnan(cf) | isinf(cf)) & ~isnan(z);
    if any(nan_mask(:))
        warning('CharacteristicFunction:NaNOutput', ...
                'cf_GL generated NaN values for valid input frequencies. Enforcing analytical limit (0) via Quant Patch.');
    end
    
    % --- Step 5: QUANT PATCH: Convert analytical collapses (NaNs) into physical limits (Zeros) ---
    cf(nan_mask) = 0;
end