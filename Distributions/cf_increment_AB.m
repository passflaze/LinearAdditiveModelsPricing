function cf_inc = cf_increment_AB(u, params, scale_factor, fwd_factor)
% CF_INCREMENT_AB  CF of the AB model increment over [t1, t2].
%
%   By additivity the t1->t2 increment of the T2-forward is
%       W = f_{T2,T2} - f_{T1,T2},   f_{T1,T2} = fwd_factor * f_{T1,T1}  (Lemma 2)
%   so its CF is the Lemma-2 rescaled ratio
%       phi_{t1->t2}(u) = phi_t2(u) / phi_t1(fwd_factor * u).
%   With fwd_factor = 1 this reduces to the plain marginal ratio phi_t2/phi_t1.
%   If scale_factor(1) == 0 (t1 = valuation date), returns phi_t2 directly.
%   NaN/Inf values from near-zero division are replaced by 0.
%
% INPUTS:
%   u            - (vector) Fourier argument
%   params       - [k; eta] model parameters passed to cf_AB
%   scale_factor - [sigma_t1*sqrt(t1), sigma_t2*sqrt(t2)] scaling factors
%   fwd_factor   - (optional, default 1) Lemma-2 rescaling B(0,T1)/B(0,T2),
%                  applied to the t1 (denominator) CF argument: phi_t1(fwd*u).
%
% OUTPUT:
%   cf_inc       - (vector) CF of the increment over [t1, t2]

    if nargin < 4 || isempty(fwd_factor)
        fwd_factor = 1;
    end

    % Special case: t1 coincides with the valuation date (pure marginal at t2)
    if scale_factor(1) == 0
        cf_inc = cf_AB(u, params, scale_factor(2));
        return;
    end

    % General case: Lemma-2 rescaled ratio  phi_t2(u) / phi_t1(fwd_factor * u)
    cf_t2  = cf_AB(u,              params, scale_factor(2));
    cf_t1  = cf_AB(fwd_factor * u, params, scale_factor(1));
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