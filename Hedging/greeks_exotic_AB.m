function g = greeks_exotic_AB(type, params_AB, scale_factor, mkt, mc, bumps)
% GREEKS_EXOTIC_AB  Delta / Gamma / Vega of an AB exotic by bump-and-revalue.
%
%   Central finite differences on price_exotic_AB. All revaluations share the
%   same random stream (mc.seed) through price_exotic_AB, so the MC noise
%   cancels in the differences.
%
%   Risk factors:
%     - Delta, Gamma : w.r.t. the underlying front future / forward F(t0,T2),
%                      bumped on mkt.forward by +/- bumps.dF (in $).
%     - Vega         : w.r.t. a PARALLEL shift of the ATM vol. Since
%                      scale_factor = sigma_t .* sqrt(yf) with sigma_t prop. to
%                      sigma_ATM, a multiplicative bump scale_factor*(1 +/- dSig)
%                      reproduces a relative sigma_ATM shift on BOTH legs.
%
% INPUTS:
%   type, params_AB, scale_factor, mkt, mc  - see price_exotic_AB
%   bumps        - (struct) .dF   absolute forward bump (e.g. 0.5 $)
%                           .dSig relative vol bump     (e.g. 1e-2)
%
% OUTPUT:
%   g - (struct) .price  base price
%                .delta  dV/dF
%                .gamma  d2V/dF2
%                .vega   dV/dsigma  (per unit of sigma_ATM; = dV/dSig / sigma)
%
%   NOTE on vega units: the bump is relative (dSig), so dV/d(relative) =
%   [V(+) - V(-)]/(2 dSig). Divide by the vol level if a per-absolute-vol
%   vega is wanted; keep the SAME convention in greeks_vanilla_AB so the
%   hedge ratios are consistent (the units cancel in build_hedge_AB).

    if nargin < 6 || isempty(bumps)
        bumps = struct('dF', 0.5, 'dSig', 1e-2);
    end

    % --- Base price -------------------------------------------------------
    V0 = price_exotic_AB(type, params_AB, scale_factor, mkt, mc);

    % --- Delta & Gamma : bump the forward ---------------------------------
    mkt_up = mkt;  mkt_up.forward = mkt.forward + bumps.dF;
    mkt_dn = mkt;  mkt_dn.forward = mkt.forward - bumps.dF;
    V_Fup = price_exotic_AB(type, params_AB, scale_factor, mkt_up, mc);
    V_Fdn = price_exotic_AB(type, params_AB, scale_factor, mkt_dn, mc);

    delta = (V_Fup - V_Fdn) / (2 * bumps.dF);
    gamma = (V_Fup - 2*V0 + V_Fdn) / (bumps.dF^2);

    % --- Vega : multiplicative bump on the (normalized) scale -------------
    V_Sup = price_exotic_AB(type, params_AB, scale_factor*(1 + bumps.dSig), mkt, mc);
    V_Sdn = price_exotic_AB(type, params_AB, scale_factor*(1 - bumps.dSig), mkt, mc);

    vega = (V_Sup - V_Sdn) / (2 * bumps.dSig);   % per relative vol bump

    % --- Pack -------------------------------------------------------------
    g = struct('price', V0, 'delta', delta, 'gamma', gamma, 'vega', vega);
end
