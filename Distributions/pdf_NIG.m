function f = pdf_NIG(x, alpha, beta, mu, delta)
% NIGPDF Calculates the Probability Density Function (PDF) of the Normal Inverse Gaussian.
%
% INPUT:
%   x     : Points at which to evaluate the PDF (scalar or vector)
%   alpha : Shape parameter (steepness/tail heaviness, alpha > |beta|)
%   beta  : Asymmetry parameter (skewness)
%   mu    : Location parameter (drift)
%   delta : Scale parameter (dispersion, delta > 0)
%
% OUTPUT:
%   f     : PDF value evaluated at x

    % 1. Calculate the recurring root term
    root_term = sqrt(delta^2 + (x - mu).^2);
    
    % 2. Constant term
    constant = (alpha * delta) / pi;
    
    % 3. Exponential term
    exponential = exp(delta * sqrt(alpha^2 - beta^2) + beta * (x - mu));
    
    % 4. Modified Bessel function of the second kind term
    bessel_term = besselk(1, alpha * root_term) ./ root_term;
    
    % 5. Assemble the final PDF
    f = constant * exponential .* bessel_term;
    
    % Protection against extreme numerical instabilities (0 * Inf)
    f(isnan(f)) = 0; 
end