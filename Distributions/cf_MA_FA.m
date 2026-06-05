function cf_inc = cf_MA_FA(u, params, scale_factor, fwd_factor)
% CF_MA_FA  CF of the Finite-Activity (FA) increment under the MA model.
%
%   By independent increments, the CF of Delta_f = f_t - f_s is the ratio of
%   the marginal CFs at t and s:
%       phi(u) = c * [(ps_plus - iu)(ps_minus + iu)] / [(pt_plus - iu)(pt_minus + iu)]
%                  * exp(i * Delta_mu * u)
%   with c = (pt_plus * pt_minus) / (ps_plus * ps_minus).
%
% INPUTS:
%   u            - (vector) Fourier argument
%   params       - [alpha; beta] MA shape parameters
%   scale_factor - [scale_s, scale_t] integrated vols at s and t
%   fwd_factor   - (optional, default 1) Lemma-2 rescaling B(0,T1)/B(0,T2);
%                  applied to the s-date CF argument: phi_s(fwd_factor * u)
% OUTPUT:
%   cf_inc       - (vector) CF of the FA increment


    % Lemma-2 forward rescaling: the t1 (s-date) marginal enters as phi_s(fwd*u),
    % i.e. u -> fwd_factor*u in every s-date factor. fwd_factor = B(0,T1)/B(0,T2);
    % default 1 (no rescaling) reproduces the plain marginal ratio phi_t/phi_s.
    if nargin < 4 || isempty(fwd_factor)
        fwd_factor = 1;
    end

    scale_factor(1) = scale_factor(1)*fwd_factor;
    
    ps_plus = params(2)/scale_factor(1); ps_minus = params(1)/scale_factor(1);
    pt_plus = params(2)/scale_factor(2); pt_minus = params(1)/scale_factor(2);
    gamma_MA = (1 / params(1)) - (1 /params(2));
    deltamu = gamma_MA * (scale_factor(2) - scale_factor(1));

    c = (pt_plus * pt_minus) / (ps_plus * ps_minus);

    numerator   = (ps_plus - 1i * u) .* (ps_minus + 1i * u);
    denominator = (pt_plus - 1i * u) .* (pt_minus + 1i * u);

    phase_shift = exp(1i * deltamu .* u);

    cf_inc = c .* (numerator ./ denominator) .* phase_shift;

    nan_mask = isnan(cf_inc) | isinf(cf_inc);
    n_bad    = nnz(nan_mask);

    if n_bad > 0
        n_nan = nnz(isnan(cf_inc));
        n_inf = nnz(isinf(cf_inc));
        warning('cf_MA_FA:NumericalInstability', ...
            '%d non-finite values in incremental CF (NaN: %d, Inf: %d) out of %d. Set to 0.', ...
            n_bad, n_nan, n_inf, numel(cf_inc));
    end

    cf_inc(nan_mask) = 0;
    
end