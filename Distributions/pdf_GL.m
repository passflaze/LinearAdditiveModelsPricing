function y = pdf_GL(params, x)
% PDF_GL  Generalized Logistic PDF, evaluated safely in log-space to prevent
%   overflow on heavy tails.
%
% INPUTS:
%   params - [alpha, beta] GL shape parameters (> 0)
%   x      - (vector) evaluation points
% OUTPUT:
%   y      - (vector) PDF values at x

    alpha    = params(1);
    beta     = params(2);
    gamma_GL = psi(beta) - psi(alpha);

    % log normalization: log(Gamma(alpha+beta) / (Gamma(alpha)*Gamma(beta)))
    log_gamma_factor = gammaln(alpha + beta) - gammaln(alpha) - gammaln(beta);

    u = x - gamma_GL;

    log_second_factor = alpha * u - (alpha + beta) * log(1 + exp(u));

    log_y = log_gamma_factor + log_second_factor;
    y = exp(log_y);


end