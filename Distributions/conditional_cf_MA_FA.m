function y = conditional_cf_MA_FA(u, params, scale_factor)
% CONDITIONAL_CF_MA_FA  Conditional-on-jump CF of the MA finite-activity increment.
%
%   Returns the CF of the MA increment over [s, t] CONDITIONAL on at least one
%   jump: (full increment CF - atom) / (1 - atom), where the no-jump atom is
%   c = (pt_plus * pt_minus) / (ps_plus * ps_minus).
%   Vanishes as |u| -> inf, so it can be safely FFT-inverted; the discrete atom
%   is handled separately by the Bernoulli split in FA_simulation.
%
% INPUTS:
%   u            - (vector) real-valued frequencies
%   params       - [alpha; beta] MA shape parameters
%   scale_factor - [scale_s; scale_t] scale factors at s (t1) and t (t2)
% OUTPUT:
%   y            - (complex vector) conditional CF evaluated at u

    ps_plus  = params(2) / scale_factor(1);
    ps_minus = params(1) / scale_factor(1);
    pt_plus  = params(2) / scale_factor(2);
    pt_minus = params(1) / scale_factor(2);

    lambda_t_s = pointmasscalculator(pt_plus, pt_minus, ps_plus, ps_minus);

    term1 = (ps_plus  - 1i.*u) ./ (pt_plus  - 1i.*u);
    term2 = (ps_minus + 1i.*u) ./ (pt_minus + 1i.*u);

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