function y = pdf_GL(params, x)
% PDF_GL_FINAL Computes the Generalized Logistic PDF safely for vectors 
% and heavy-tails using log-space computation to prevent overflow.

    % 1. Location parameter (Difference of Digamma functions)
    alpha = params(1); beta = params(2);
    gamma_GL = psi(beta) - psi(alpha); 
    
    % 2. Logarithm of the normalizing Gamma/Beta factor
    % Written as log(gamma(a+b)/(gamma(a)*gamma(b))) using gammaln for stability
    log_gamma_factor = gammaln(alpha + beta) - gammaln(alpha) - gammaln(beta);

    u = x - gamma_GL;
    
    % Body of the distribution (Standard safe evaluation)
    log_second_factor = alpha * u - ...
        (alpha + beta) * log(1 + exp(u));


    % 4. Combine in log-space and exponentiate at the very end
    log_y = log_gamma_factor + log_second_factor;
    y = exp(log_y);


end