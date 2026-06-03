function g = greeks_vanilla_AB(instrument, K, maturity_index, params_AB, scale_factors, mkt, params_hedge, mc, bumps)
% GREEKS_VANILLA_AB  Price + Delta / Gamma / Vega of a SINGLE hedging vanilla.
%
%   Building block of a multi-instrument vanilla hedge: a plain-vanilla call,
%   a plain-vanilla put or the underlying future, at an explicit strike K and
%   maturity leg. Greeks use the SAME bump definitions as greeks_exotic_AB so
%   that the ratios in build_hedge are consistent (all quantities are t0 PV).
%
%       future : price = F,  delta = 1,  gamma = 0,  vega = 0
%       call   : AB analytic price (Lewis-FFT, discounted to t0), FD Greeks
%       put    : via put-call parity  P = C + B * (K - F)
%
% INPUTS:
%   instrument     - (string) 'call' | 'put' | 'future'
%   K              - (scalar) strike (ignored for 'future')
%   maturity_index - (scalar) maturity leg index (e.g. 4 for T2)
%   params_AB      - (2x1)    [k; eta]
%   scale_factors  - (vector) sigma_t .* sqrt(yf) per maturity
%   mkt            - (struct) Market data
%   params_hedge   - (struct) Hedging parameters (forward, ...)
%   mc             - (struct) Monte Carlo/FFT parameters
%   bumps          - (struct) .dF .dSig   (same as greeks_exotic_AB)
%
% OUTPUT:
%   g - (struct) .kind .K .mat .price .delta .gamma .vega

    if nargin < 9 || isempty(bumps)
        bumps = struct('dF', 0.5, 'dSig', 1e-2);
    end

    B  = mkt.discount_factor(maturity_index); % discount to maturity
    F  = mkt.forward(maturity_index);
    scale_factor = scale_factors(maturity_index);

    switch lower(instrument)
        case 'future'
            g = struct('kind', 'future', 'K', NaN, 'mat', maturity_index, ...
                       'price', F, 'delta', 1, 'gamma', 0, 'vega', 0);
            return;
        case 'call'
            isCall = true;
        case 'put'
            isCall = false;
        otherwise
            error('greeks_vanilla_AB:badInstrument', 'Unknown instrument "%s".', instrument);
    end

    % --- Base + bumped CALL prices (central FD), discounted to t0 ----------
    % lewis_FFT_call returns the UNDISCOUNTED (forward) call value, so we
    % multiply by B to get a t0 PV in the SAME measure as the exotic prices.
    % The 7th argument (=1) is the doubleshift flag; fwd_factor (9th) = 1.
    C0    = B * lewis_FFT_call(@cf_AB, mc.M, mc.dz, params_AB, scale_factor, K - F,            1, 'AB');
    C_Fup = B * lewis_FFT_call(@cf_AB, mc.M, mc.dz, params_AB, scale_factor, K - (F+bumps.dF), 1, 'AB');
    C_Fdn = B * lewis_FFT_call(@cf_AB, mc.M, mc.dz, params_AB, scale_factor, K - (F-bumps.dF), 1, 'AB');

    % --- Vega (Call Vega = Put Vega at the same strike) -------------------
    vega = compute_vega_AB('vanilla', mkt, params_hedge, mc, bumps.dSig, maturity_index, K);

    % --- Put-Call Parity Adjustments (all quantities are t0 PV) -----------
    if isCall
        V0    = C0;
        V_Fup = C_Fup;
        V_Fdn = C_Fdn;
    else
        % Put-Call Parity in t0 PV: P = C + B(t0,T) * (K - F)
        V0    = C0    + B * (K - F);
        V_Fup = C_Fup + B * (K - (F + bumps.dF));
        V_Fdn = C_Fdn + B * (K - (F - bumps.dF));
    end

    % --- Final Greeks Assembly --------------------------------------------
    g = struct( ...
        'kind',  lower(instrument), ...
        'K',     K, ...
        'mat',   maturity_index, ...
        'price', V0, ...
        'delta', (V_Fup - V_Fdn) / (2 * bumps.dF), ...
        'gamma', (V_Fup - 2 * V0 + V_Fdn) / (bumps.dF^2), ...
        'vega',  vega);
end
