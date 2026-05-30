function y = cf_FA_MA(u, pt_plus, pt_minus, ps_plus, ps_minus)
%CF_FA_MA Computes the conditional characteristic function of a single jump.
%   Y = CF_FA_MA(U, PT_PLUS, PT_MINUS, PS_PLUS, PS_MINUS) calculates the 
%   characteristic function of individual jumps, phi_J(u), conditional on 
%   at least one jump occurring within the Minimal Additive (MA) framework.
%
%   This conditional characteristic function is crucial for finite-activity 
%   Lévy-driven formulations because it naturally vanishes at infinity 
%   (lim_{|u|->inf} y = 0), ensuring the numerical stability and 
%   applicability of FFT-based inversion methods.
%
%   The function is fully vectorized and supports both scalar and array 
%   inputs for the frequency variable U.
%
%   Inputs:
%       u        - Scalar or array of real-valued frequencies (grid points)
%       pt_plus  - Right-tail decay parameter at terminal time t (pt_plus > 0)
%       pt_minus - Left-tail decay parameter at terminal time t (pt_minus > 0)
%       ps_plus  - Right-tail decay parameter at initial time s (ps_plus > 0)
%       ps_minus - Left-tail decay parameter at initial time s (ps_minus > 0)
%
%   Outputs:
%       y        - Complex characteristic function values evaluated at each U
%

    % 1. Compute the integrated intensity Lambda over [s, t]
    lambda_t_s = pointmasscalculator(pt_plus, pt_minus, ps_plus, ps_minus);

    % 2. Vectorized Rational Fraction Terms
    term1 = (ps_plus - 1i.*u) ./ (pt_plus - 1i.*u);
    term2 = (ps_minus + 1i.*u) ./ (pt_minus + 1i.*u);

    % 3. Conditional CF Calculation
    phi = (1 / lambda_t_s) .* log(term1 .* term2);

    num = exp(lambda_t_s.*phi)-1;
    den = exp(lambda_t_s)-1;
    y = num./den;

    nan_mask = isnan(y);
    inf_mask = isinf(y);
    
    if any(nan_mask) || any(inf_mask)
        num_nan = sum(nan_mask(:));
        num_inf = sum(inf_mask(:));
        
        fprintf('\n[WARNING] Numerical instability detected in %s:\n', mfilename);
        if num_nan > 0
            fprintf(' -> Found %d NaN element(s) in the output.\n', num_nan);
        end
        if num_inf > 0
            fprintf(' -> Found %d Inf element(s) in the output.\n', num_inf);
        end
        fprintf('          Check for identical boundary parameters (e.g., pt_plus == ps_plus) or extreme U grid ranges.\n\n');
    end
    
end