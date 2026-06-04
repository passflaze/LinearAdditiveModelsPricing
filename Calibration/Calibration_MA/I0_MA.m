function I0 = I0_MA(params)
% I0_MA Normalization constant I0 of the Minimal Additive (MA) model.
%
% I0 ties the market ATM volatility to the internal normalized scale:
%       sigma_t = sigma_ATM / I0,
% with I0 = sqrt(2*pi) * E[z^+] for z the standardized MA marginal. The
% closed-form expression is piecewise in the sign of the centering drift
% gamma_MA = 1/alpha - 1/beta (alpha and beta are derived from params).
%
% INPUTS:
%   params : [alpha; beta] MA shape parameters (both strictly non-zero)
%
% OUTPUT:
%   I0     : analytical normalization constant (scalar), a pure number
%            depending only on (alpha, beta)


    alpha = params(1); beta = params(2);
    if alpha == 0 || beta == 0
        error('I0_MA:ZeroParameter', 'alpha and beta must be strictly non-zero.');
    end

    gamma_MA = (1/alpha) - (1/beta);
    C        = 1 / ((1/beta) + (1/alpha));

    % Piecewise on the sign of the centering drift gamma_MA.
    if gamma_MA > 0
        I0 = sqrt(2*pi) * (C / (alpha^2)) * exp(-gamma_MA * alpha);
    else
        I0 = sqrt(2*pi) * (C / (beta^2)) * exp(gamma_MA * beta);
    end
end