function [price] = call_ATM_vanilla(params, scale_factor_T, discount_factor_T, model)
% CALL_ATM_VANILLA  Closed-form price of an ATM-forward European vanilla
%                   call under the GL or AB additive model.
%
%   Under martingality of the forward (zero-mean increment), an ATM-forward
%   call has the dimensionally homogeneous closed form:
%
%       C_ATM(T) = B(0,T) * sigma_t(T) * sqrt(T) * normalized_ATM_price
%
%   where:
%     - B(0,T)                  : discount factor at T
%     - sigma_t(T) * sqrt(T)    : NORMALIZED scale (= input scale_factor_T)
%     - normalized_ATM_price    : E[z^+], the expected positive part of the
%                                 standardized model variable z. This is a
%                                 pure number depending only on the model's
%                                 shape parameters.
%
%   We do NOT need a separate utility for the normalized ATM price: it is
%   exactly  I_0 / sqrt(2*pi), where I_0 is the same constant already used
%   to define the normalized scale  sigma_t = sigma_ATM / I_0. So we just
%   reuse the existing I0_GL / I0_AB.
%
%   No Monte Carlo, no FFT: this is the analytic reference value used as
%   sanity check against MC-based pricers (e.g. the K1 = 0 collapse of
%   CoC_pricing_MC must reproduce this value).
%
% INPUTS:
%   params              - (1x2) model parameters, same convention as in
%                         CoC_pricing_MC:
%                           * 'GL': [alpha, beta]
%                           * 'AB': [eta,   kappa]
%   scale_factor_T      - (scalar) sigma_t(T) * sqrt(T) at the call maturity.
%                         Must be the SAME normalized scale used inside
%                         CoC_pricing_MC and Chooser_pricing_MC.
%   discount_factor_T   - (scalar) B(0,T) at the call maturity.
%   model               - (string) 'GL' or 'AB'.
%
% OUTPUTS:
%   price    - (scalar) ATM-forward vanilla call price at t0.
%   details  - (struct) audit info:
%                .normalized_ATM_price   E[z^+], pure number
%                .scale_factor_T         echo of input
%                .discount_factor_T      echo of input
%                .model                  echo of input


    % --- Normalized ATM price  E[z^+] = I_0 / sqrt(2*pi) -----------------
    switch model
        case 'GL'
            I0_value = I0_GL(params);
        case 'AB'
            % I0_AB convention in this codebase: I0_AB(x0, params)
            % evaluated at the moneyness x0 = 0 (ATM normalized).
            I0_value = I0_AB(0, params);
        otherwise
            error('call_ATM_vanilla:badModel', ...
                'Unknown model "%s". Use ''GL'' or ''AB''.', model);
    end

    normalized_ATM_price = I0_value / sqrt(2*pi);

    % --- Dimensional price -----------------------------------------------
    price = discount_factor_T * scale_factor_T * normalized_ATM_price;

end