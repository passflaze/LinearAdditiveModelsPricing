function I0 = I0_GL(params)
% I0_GL  Normalization constant I_0 of the GL model.
%
%   I_0 ties the market ATM volatility sigma_ATM to the internal normalized
%   scale sigma_t used by the pricers:    sigma_t = sigma_ATM / I_0.
%
%   Definition (same as in the main script):
%       I_0 = sqrt(2*pi) * integral_0^inf  x * pdf_GL(params, x) dx
%
%   Equivalently  I_0 = sqrt(2*pi) * E[z^+]  with z ~ standardized GL.
%   The sqrt(2*pi) factor comes from the Bachelier ATM-call formula
%   C_ATM = S * sqrt(T/(2*pi)) * sigma_ATM, which fixes the calibration
%   of sigma_t to market quotes.
%
% INPUT:
%   params - (1x2) GL shape parameters [alpha, beta]
%
% OUTPUT:
%   I0     - (scalar) normalization constant, pure number depending only
%            on (alpha, beta).

    integrand_positive_part = @(x) pdf_GL(params, x) .* x;
    I0 = sqrt(2*pi) * quadgk(integrand_positive_part, 0, inf);
end