function y = conditional_cf_MA_FA(u, params, scale_factor)
%CONDITIONAL_CF_MA_FA Computes the conditional characteristic function of a single jump.
%   Y = CONDITIONAL_CF_MA_FA(U, PT_PLUS, PT_MINUS, PS_PLUS, PS_MINUS) calculates the 
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

    ps_plus = params(2)/scale_factor(1); ps_minus = params(1)/scale_factor(1);
    pt_plus = params(2)/scale_factor(2); pt_minus = params(1)/scale_factor(2);

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
    
    n_bad = nnz(nan_mask) + nnz(inf_mask);
    if n_bad > 0
        % Use warning() (not fprintf): silenceable and non-spamming inside MC
        % loops. Identical boundary params (e.g. pt_plus == ps_plus) or an
        % extreme u-grid are the usual cause.
        warning('conditional_cf_MA_FA:NumericalInstability', ...
            ['%d non-finite value(s) (NaN: %d, Inf: %d) out of %d. ', ...
             'Check for identical boundary parameters or extreme u-grid ranges.'], ...
            n_bad, nnz(nan_mask), nnz(inf_mask), numel(y));
    end

end