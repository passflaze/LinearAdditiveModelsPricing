function cf = cf_GL(u, params, scale_factor)
% CF_GL  Characteristic function of the Generalized Laplace distribution.
%
%   phi(u) = [Gamma(alpha + i*u) * Gamma(beta - i*u) / Gamma(alpha+beta)]
%            / B(alpha, beta) * exp(i * gamma_GL * u)
%
%   where gamma_GL = psi(beta) - psi(alpha) is the centering drift.
%   NaN/Inf values are replaced by 0 (tail regularisation).
%
% INPUTS:
%   u            - (vector) Fourier argument
%   params       - [alpha, beta] shape parameters (must be real and > 0)
%   scale_factor - (scalar) time/volatility scaling applied to u
%
% OUTPUT:
%   cf           - (vector) CF evaluated at u

    alpha = params(1);
    beta  = params(2);

    % Validate: shape parameters must be real
    if any(imag(alpha) ~= 0) || any(imag(beta) ~= 0)
        error('cf_GL:InvalidParameters', ...
            'Shape parameters alpha and beta must be real-valued.');
    end

    % Centering drift (ensures zero mean)
    gamma_GL = psi(beta) - psi(alpha);

    % Characteristic function via Gamma function ratio
    num = complex_gamma(alpha + 1i .* u * scale_factor) .* ...
          complex_gamma(beta  - 1i .* u * scale_factor);
    den = complex_gamma(alpha + beta) .* exp(betaln(alpha, beta));

    cf  = (num ./ den) .* exp(1i .* gamma_GL .* u * scale_factor);

    % Numerical stability: replace NaN/Inf from tail collapse with 0
    nan_mask = (isnan(cf) | isinf(cf)) & ~isnan(u);

    if any(nan_mask(:))
        n_bad = nnz(nan_mask);
        warning('cf_GL:NumericalInstability', ...
            '%d non-finite values in CF (out of %d). Set to 0.', ...
            n_bad, numel(cf));
    end

    cf(nan_mask) = 0;

end