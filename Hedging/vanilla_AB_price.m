function price = vanilla_AB_price(F, K, scale_factor_T, B, params_AB, isCall)
% Dollar price of a single AB vanilla, in the SAME scale convention used by
% the exotic pricers (scale_factor_T = sigma_t*sqrt(T) = sigma_ATM*sqrt(T)/I0).
%
%   sigma_ATM*sqrt(T) = scale_factor_T * I0
%   chi               = (K - F) / (sigma_ATM*sqrt(T))
%   C(K)              = B * sigma_ATM*sqrt(T) * G(chi)         (Baviera-Massaria Eq.20)
%
% TODO(Persona A): sanity-check this against call_ATM_vanilla at K = F:
%   vanilla_AB_price(F,F,s,B,p,true) must equal call_ATM_vanilla(p,s,B,'AB').
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