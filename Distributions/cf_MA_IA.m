function cf = cf_MA_IA(u, params, scale_factors)
% CF_MA_IA  Characteristic function of the MA marginal (base process at time t).
%
%   Evaluates the exact marginal CF of the Minimal Additive process. Tail decays
%   are derived from the largest scale factor: pt_plus = beta/scale,
%   pt_minus = alpha/scale; gamma_t centres the process (zero mean).
%
% INPUTS:
%   u             - (vector) Fourier argument
%   params        - [alpha; beta] MA shape parameters
%   scale_factors - scale factor(s); max(scale_factors) = sigma_t * sqrt(t) is used
% OUTPUT:
%   cf            - (complex array) CF evaluated at u

    scale_factor = max(scale_factors);
    pt_plus  = params(2) / scale_factor;
    pt_minus = params(1) / scale_factor;

    gamma_t = (1 ./ pt_minus) - (1 ./ pt_plus);

    % Log-CF: -log(1 - iu*(1/pt_plus - 1/pt_minus) + u^2/(pt_plus*pt_minus))
    term_linear = 1i .* u .* ((1 ./ pt_plus) - (1 ./ pt_minus));
    term_quad   = (u.^2) ./ (pt_plus .* pt_minus);
    poly_term   = 1 - term_linear + term_quad;

    lcf = -log(poly_term) + 1i .* gamma_t .* u;

    cf = exp(lcf);

end