function cf = cf_GL_cal(z, alpha, beta)
% CF_GL_CAL  Characteristic function of the Generalized Logistic (unit scale).
%   phi(z) = [Gamma(alpha+i z) Gamma(beta-i z) / (Gamma(alpha) Gamma(beta))]
%            * exp(i*gamma*z),  gamma = psi(beta) - psi(alpha) (zero-mean drift).
%
% INPUTS:
%   z          : frequencies (real or complex)
%   alpha, beta: shape parameters (> 0, real). alpha = left rate, beta = right.
% OUTPUT:
%   cf         : characteristic function at z (non-finite tails set to 0).

    if any(imag(alpha) ~= 0) || any(imag(beta) ~= 0)
        error('cf_GL_cal:ComplexParameters', 'alpha and beta must be real.');
    end

    gamma_GL = psi(beta) - psi(alpha);

    num = (complex_gamma(alpha + 1i.*z) .* complex_gamma(beta - 1i.*z)) ./ complex_gamma(alpha + beta);
    den = exp(betaln(alpha, beta));

    cf = (num ./ den) .* exp(1i .* gamma_GL .* z);

    % Set non-finite tail values to their analytical limit (0).
    cf((isnan(cf) | isinf(cf)) & ~isnan(z)) = 0;
end