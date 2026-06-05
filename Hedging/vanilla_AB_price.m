function price = vanilla_AB_price(F, K, scale_factor_T, B, params_AB, isCall)
% VANILLA_AB_PRICE  Dollar price of a single AB vanilla option or put.
%   Uses the SAME scale convention as the exotic pricers:
%       scale_factor_T = sigma_t * sqrt(T) = sigma_ATM * sqrt(T) / I0
%       chi = (K - F) / (sigma_ATM * sqrt(T))
%       C(K) = B * sigma_ATM*sqrt(T) * G(chi)    (Baviera-Massaria Eq. 20)
%
% INPUTS:
%   F              - (scalar) forward price F(t0, T)
%   K              - (scalar) strike
%   scale_factor_T - (scalar) sigma_t * sqrt(T) at maturity T
%   B              - (scalar) discount factor B(t0, T)
%   params_AB      - (2x1)    [k; eta]
%   isCall         - (logical) true for call, false for put (put-call parity)
% OUTPUT:
%   price          - (scalar) discounted option price at t0
    k   = params_AB(1);
    eta = params_AB(2);
    I0  = I0_AB(0, params_AB);

    sigATM_sqrtT = scale_factor_T * I0;
    chi          = (K - F) / sigATM_sqrtT;
    G            = lewis_FFT_call(chi, k, eta, I0);
    call         = B * sigATM_sqrtT * G;

    if isCall
        price = call;
    else
        price = call - (F - K) * B;   % put-call parity in the forward measure
    end
end