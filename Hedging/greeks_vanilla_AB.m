function g = greeks_vanilla_AB(instrument, params_AB, scale_factor_T, mkt, bumps)
% GREEKS_VANILLA_AB  Price + Delta / Gamma / Vega of a hedging instrument.
%
%   Building blocks of the hedge: plain-vanilla call, plain-vanilla put and
%   the underlying future, all on maturity T2. Greeks use the SAME bump
%   definitions as greeks_exotic_AB so that the ratios in build_hedge_AB are
%   consistent.
%
%       future : price = F,  delta = 1,  gamma = 0,  vega = 0
%       call   : AB analytic price (Lewis-FFT, call_AB_FFT), FD Greeks
%       put    : via put-call parity  P = C - (F - K) * B
%
% INPUTS:
%   instrument     - (string) 'call' | 'put' | 'future'
%   params_AB      - (2x1) [k; eta]
%   scale_factor_T - (scalar) sigma_t(T2)*sqrt(T2)  (same scale as exotic, leg 2)
%   mkt            - (struct) .forward F(t0,T2)
%                             .Kc / .Kp strikes of the hedge call / put
%                             .df    [B(t0,T1), B(t0,T2)]  (uses B(t0,T2))
%   bumps          - (struct) .dF .dSig   (same as greeks_exotic_AB)
%
% OUTPUT:
%   g - (struct) .price .delta .gamma .vega

    if nargin < 5 || isempty(bumps)
        bumps = struct('dF', 0.5, 'dSig', 1e-2);
    end

    B  = mkt.df(2);          % discount to T2
    F  = mkt.forward;

    switch lower(instrument)
        case 'future'
            g = struct('price', F, 'delta', 1, 'gamma', 0, 'vega', 0);
            return;

        case 'call'
            K = mkt.Kc;  isCall = true;
        case 'put'
            K = mkt.Kp;  isCall = false;
        otherwise
            error('greeks_vanilla_AB:badInstrument', ...
                'Unknown instrument "%s".', instrument);
    end

    % --- Base + bumped prices (central FD), CRN not needed (analytic) -----
    V0    = vanilla_AB_price(F,             K, scale_factor_T,              B, params_AB, isCall);
    V_Fup = vanilla_AB_price(F + bumps.dF,  K, scale_factor_T,              B, params_AB, isCall);
    V_Fdn = vanilla_AB_price(F - bumps.dF,  K, scale_factor_T,              B, params_AB, isCall);
    V_Sup = vanilla_AB_price(F,             K, scale_factor_T*(1+bumps.dSig), B, params_AB, isCall);
    V_Sdn = vanilla_AB_price(F,             K, scale_factor_T*(1-bumps.dSig), B, params_AB, isCall);

    g = struct( ...
        'price', V0, ...
        'delta', (V_Fup - V_Fdn) / (2*bumps.dF), ...
        'gamma', (V_Fup - 2*V0 + V_Fdn) / (bumps.dF^2), ...
        'vega',  (V_Sup - V_Sdn) / (2*bumps.dSig));
end

% =========================================================================
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
    G            = call_AB_FFT(chi, k, eta, I0);
    call         = B * sigATM_sqrtT * G;

    if isCall
        price = call;
    else
        price = call - (F - K) * B;   % put-call parity in the forward measure
    end
end
