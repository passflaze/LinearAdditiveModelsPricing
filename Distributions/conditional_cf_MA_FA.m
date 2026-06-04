function y = conditional_cf_MA_FA(u, params, scale_factor)
%CONDITIONAL_CF_MA_FA Conditional-on-jump CF of the MA finite-activity increment.
%   Y = CONDITIONAL_CF_MA_FA(U, PARAMS, SCALE_FACTOR) returns the CF of the MA
%   increment over [s, t] CONDITIONAL on at least one jump occurring, i.e.
%   (full increment CF - atom)/(1 - atom), where the atom is the no-jump point
%   mass c = (pt_plus*pt_minus)/(ps_plus*ps_minus). The tail decays are derived
%   internally from the parameters and the two scale factors.
%
%   This conditional CF vanishes at infinity (lim_{|u|->inf} y = 0), so it can
%   be safely FFT-inverted, while the discrete atom is handled separately by the
%   Bernoulli split in FA_simulation (cf. project hint 3.c.i and [5]).
%
%   The function is fully vectorized in the frequency variable U.
%
%   Inputs:
%       u            - real-valued frequencies (scalar or array)
%       params       - [alpha; beta] MA shape parameters
%       scale_factor - [scale_s; scale_t] scale factors at s (t1) and t (t2)
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