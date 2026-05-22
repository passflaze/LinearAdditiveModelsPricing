function pdf_GL = pdf_GL(alpha, beta, x)
%PDF_GL Summary of this function goes here
%   Detailed explanation goes here


% Translation parameter (using MATLAB's built-in digamma function 'psi')
gamma_GL = psi(beta) - psi(alpha); % to check !!!!

% Normalization constant
C_GL = gamma(alpha + beta) / (gamma(alpha) * gamma(beta));

% PDF definition (fully vectorized)
pdf_GL = C_GL .* (exp(alpha .* (x - gamma_GL)) ./ (1 + exp(x - gamma_GL)).^(alpha + beta));
end