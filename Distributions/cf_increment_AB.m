function cf_inc = cf_increment_AB(u, params, scale_factor)
% CF_INCREMENT_AB  CF of the AB model increment over [t1, t2].
%
%   Uses the additivity property: phi_{t1->t2}(u) = phi_t2(u) / phi_t1(u).
%   If scale_factor(1) == 0 (t1 = valuation date), returns phi_t2 directly.
%   NaN/Inf values from near-zero division are replaced by 0.
%
% INPUTS:
%   u            - (vector) Fourier argument
%   params       - [eta, kappa] model parameters passed to cf_AB
%   scale_factor - [sigma_t1*sqrt(t1), sigma_t2*sqrt(t2)] scaling factors
%
% OUTPUT:
%   cf_inc       - (vector) CF of the increment over [t1, t2]

    % Special case: t1 coincides with the valuation date
    if scale_factor(1) == 0
        cf_inc = cf_AB(u, params, scale_factor(2));
        return;
    end

    % General case: ratio of the two cumulative CFs
    cf_t2  = cf_AB(u, params, scale_factor(2));  
    cf_t1  = cf_AB(u, params, scale_factor(1));   
    cf_inc = cf_t2 ./ cf_t1;

    % Numerical stability: replace NaN/Inf from near-zero division with 0
    nan_mask = isnan(cf_inc) | isinf(cf_inc);
    n_bad    = nnz(nan_mask);

    if n_bad > 0
        n_nan = nnz(isnan(cf_inc));
        n_inf = nnz(isinf(cf_inc));
        warning('cf_increment_AB:NumericalInstability', ...
            '%d non-finite values in incremental CF (NaN: %d, Inf: %d) out of %d. Set to 0.', ...
            n_bad, n_nan, n_inf, numel(cf_inc));
    end

    cf_inc(nan_mask) = 0;

end