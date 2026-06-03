function g = greeks_exotic_AB(type, params_AB, scale_factor, market, params_hedge, mc, bumps)
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
        bumps = struct('dF', 1e-4, 'dSig', 1e-4);
    end

    % --- Base price -------------------------------------------------------
    V0 = price_exotic_AB(type, params_AB, scale_factor, market,params_hedge, mc);

    % --- Delta & Gamma : bump the forward ---------------------------------
    params_up = params_hedge;  params_up.forward = params_hedge.forward + bumps.dF;
    params_down = params_hedge;  params_down.forward = params_hedge.forward - bumps.dF;
    V_Fup = price_exotic_AB(type, params_AB, scale_factor, market,params_up, mc);
    V_Fdn = price_exotic_AB(type, params_AB, scale_factor, market,params_down, mc);

    delta = (V_Fup - V_Fdn) / (2 * bumps.dF);
    gamma = (V_Fup - 2*V0 + V_Fdn) / (bumps.dF^2);

    
    vega = compute_vega_AB(type, market,params_hedge, mc, bumps.dSig);

    % --- Pack -------------------------------------------------------------
    g = struct('price', V0, 'delta', delta, 'gamma', gamma, 'vega', vega);
end
