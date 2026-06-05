function I0 = I0_GL(params)
% I0_GL  Normalization constant I0 = sqrt(2*pi)*E[z^+] of the GL model, with z
%        the standardized GL marginal. Ties market ATM vol to the internal
%        scale: sigma_t = sigma_ATM / I0.  (sqrt(2*pi) from the Bachelier ATM
%        call.)
%
% INPUTS:
%   params : [alpha; beta] GL shape parameters (> 0)
% OUTPUT:
%   I0     : normalization constant (scalar) depending only on (alpha, beta)

    integrand_positive_part = @(x) pdf_GL(params, x) .* x;
    I0 = sqrt(2*pi) * quadgk(integrand_positive_part, 0, inf);
end