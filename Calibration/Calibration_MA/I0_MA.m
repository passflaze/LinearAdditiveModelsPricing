function I0 = I0_MA(alpha, beta)
% I0_MA Computes the normalization constant (I0) for the MA model.
%
% This function calculates the analytical structural component I0 used to
% enforce internal consistency or martingale restrictions within the MA 
% (Moving Average / Modified Analytical) pricing framework. The calculation
% switches analytically based on the sign of the threshold parameter gamma_MA.
%
% INPUTS:
%   gamma_MA : Threshold log-moneyness parameter (scalar)
%   C        : Model scaling/intensity coefficient (must be non-zero)
%   alpha    : Left tail decay parameter (must be non-zero)
%   beta     : Right tail decay parameter (must be non-zero)
%
% OUTPUT:
%   I0       : Analytical normalization factor value (scalar)

    % =========================================================================
    % STEP 1: INPUT INTEGRITY CHECKS
    % =========================================================================
    if alpha == 0 || beta == 0
        error('I0_MA:ZeroParameter', ...
              'Parameters alpha and beta must be strictly non-zero to prevent division by zero.');
    end

    gamma_MA = (1/alpha) - (1/beta);
    C     = 1 / ((1/beta) + (1/alpha)); 
    % =========================================================================
    % STEP 2: PIECEWISE ANALYTICAL EVALUATION
    % =========================================================================
    % The structural integral changes its convergence form depending on whether
    % the log-moneyness threshold gamma_MA is positive or negative.
    if gamma_MA > 0
        % Evaluation for positive regime
        I0 = sqrt(2*pi) * (C / (alpha^2)) * exp(-gamma_MA * alpha);
    else
        % Evaluation for negative/zero regime
        I0 = sqrt(2*pi) * (C / (beta^2)) * exp(gamma_MA * beta);
    end
end