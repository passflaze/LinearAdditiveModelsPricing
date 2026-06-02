function price = vanilla_AB_price(F, K, scale_factor_T, B, params_AB, isCall)
% VANILLA_AB_PRICE  Dollar price of a single AB vanilla (call or put).
%
%   Same scale convention used by the exotic pricers:
%       scale_factor_T = sigma_t*sqrt(T) = sigma_ATM*sqrt(T)/I0
%       sigma_ATM*sqrt(T) = scale_factor_T * I0
%       chi  = (K - F) / (sigma_ATM*sqrt(T))
%       C(K) = B * sigma_ATM*sqrt(T) * G(chi)        (Baviera-Massaria Eq.20)
%       P(K) = C(K) - (F - K) * B                    (put-call parity, fwd measure)
%
%   ATM sanity: vanilla_AB_price(F,F,s,B,p,true) == call_ATM_vanilla(p,s,B,'AB').
%
% INPUTS:
%   F, K           - forward and strike (scalars)
%   scale_factor_T - sigma_t(T)*sqrt(T) at the option maturity
%   B              - discount factor B(0,T)
%   params_AB      - [k; eta]
%   isCall         - true for a call, false for a put
%
% OUTPUT:
%   price          - dollar price at t0.
    k   = params_AB(1);
    eta = params_AB(2);
    I0  = I0_AB(0, params_AB);

    sigATM_sqrtT = scale_factor_T * I0;
    chi          = (K - F) / sigATM_sqrtT;
    G            = call_AB_FFT(chi, k, eta, I0);
    call         = B * sigATM_sqrtT * G;

    if isCall
        price = call;
    else
        price = call - (F - K) * B;
    end
end
