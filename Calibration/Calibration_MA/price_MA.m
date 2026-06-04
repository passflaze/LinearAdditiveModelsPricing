function CallPrice = price_MA(params, discount_factor, yearfraction, sigma_ATM, moneyness_modified)
%PRICE_MA  European call surface under the Minimal Additive model (closed form).
%   C(K,t) = B(t)*(sigma_ATM/I0)*sqrt(t) * G(I0*chi), with G the normalized MA
%   call, piecewise around the threshold X = I0*chi - gamma_MA.
%
% INPUTS:
%   params             : [alpha, beta]  (both > 0)
%   discount_factor, yearfraction, sigma_ATM : (M x 1) per maturity
%   moneyness_modified : (M x K) chi = (K-F)/(sigma_ATM*sqrt(t))
% OUTPUT:
%   CallPrice          : (M x K) call prices

    alpha = params(1);
    beta  = params(2);

    gamma_MA = (1/alpha) - (1/beta);
    C        = 1 / ((1/beta) + (1/alpha));
    I0       = I0_MA(params);

    maturity_multiplier = discount_factor .* (sigma_ATM / I0) .* sqrt(yearfraction);   % (M x 1)
    X = I0 * moneyness_modified - gamma_MA;                                            % (M x K)

    core_value = zeros(size(moneyness_modified));
    mask_less    = (X < 0);          % ITM side (I0*chi < gamma_MA)
    mask_greater = ~mask_less;       % OTM side
    core_value(mask_less)    = (C / (alpha^2)) * exp(alpha * X(mask_less)) - moneyness_modified(mask_less)*I0;
    core_value(mask_greater) = (C / (beta^2))  * exp(-beta * X(mask_greater));

    CallPrice = maturity_multiplier .* core_value;

end