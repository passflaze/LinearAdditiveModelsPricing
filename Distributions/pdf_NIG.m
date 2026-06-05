function f = pdf_NIG(x, alpha, beta, mu, delta)
% PDF_NIG  Probability density function of the Normal Inverse Gaussian distribution.
%
% INPUTS:
%   x     - (vector) evaluation points
%   alpha - (scalar) shape / tail-heaviness parameter (alpha > |beta|)
%   beta  - (scalar) asymmetry parameter (skewness)
%   mu    - (scalar) location parameter
%   delta - (scalar) scale parameter (> 0)
% OUTPUT:
%   f     - (vector) PDF values at x

    root_term   = sqrt(delta^2 + (x - mu).^2);
    constant    = (alpha * delta) / pi;
    exponential = exp(delta * sqrt(alpha^2 - beta^2) + beta * (x - mu));
    bessel_term = besselk(1, alpha * root_term) ./ root_term;

    f = constant * exponential .* bessel_term;

    % Guard against 0 * Inf at extreme tails
    f(isnan(f)) = 0;
end