function I0 = I0_GL(alpha, beta)
% I0_GL Computes the scaled partial expectation of a Gauss-Laplace distribution.
%
%   I0 = I0_GL(ALPHA, BETA) calculates the integral of x * pdf_GL(x) 
%   over the domain [0, +inf], scaled by the constant sqrt(2*pi).
%
%   Inputs:
%       alpha - The alpha parameter of the Gauss-Laplace distribution
%       beta  - The beta parameter of the Gauss-Laplace distribution
%
%   Outputs:
%       I0    - The computed numerical integral value
%
%   Dependencies:
%       pdf_GL - Must be a valid function in your MATLAB path or defined below.

    % =========================================================================
    % INTEGRAND DEFINITION
    % =========================================================================
    % Define the integrand for the partial moment (Expected value over [0, inf])
    % We use element-wise multiplication (.*) to ensure vectorization safety 
    % within the quadgk algorithm.
    integrand_mean = @(x) pdf_GL(alpha, beta, x) .* x;

    % =========================================================================
    % NUMERICAL INTEGRATION
    % =========================================================================
    % Compute the integral using Gauss-Kronrod quadrature (quadgk)
    % which is highly efficient for infinite domains (0 to inf).
    I0 = sqrt(2*pi) * quadgk(integrand_mean, 0, inf);

end