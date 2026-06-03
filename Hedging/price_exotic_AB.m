function V = price_exotic_AB(type, params_AB, scale_factor, mkt, params_hedge, mc)
% PRICE_EXOTIC_AB  Uniform repricing wrapper for the AB exotics.
%
%   Single entry point that dispatches to the existing MC pricers
%   (CoC_pricing_MC / PoP_pricing_MC / Chooser_pricing_MC) so that the
%   Greek engine (greeks_exotic_AB) can bump-and-revalue through ONE call.
%
%   COMMON RANDOM NUMBERS: rng(mc.seed) is reset before every MC pricing,
%   so two evaluations that differ only by a small bump use the SAME random
%   stream. This is mandatory for clean finite-difference Greeks (especially
%   gamma, a second difference dominated by MC noise otherwise).
%
% INPUTS:
%   type         - (string) 'CoC' | 'PoP' | 'Chooser'
%   params_AB    - (2x1)    [k; eta]
%   scale_factor - (1x2)    [scale_t1, scale_t2]  (= sigma_t .* sqrt(yf))
%   params_hedge          - (struct) .forward  F(t0,T2)
%                           .K1       compound strike (ignored by Chooser)
%                           .K2       inner strike
%                           .df       [B(t0,T1), B(t0,T2)]
%   mc           - (struct) .N_sim .M .dz .N_grid .seed
%
% OUTPUT:
%   V            - (scalar) exotic fair value at t0 under the AB model.

    model = 'AB';

    discount_factors = [mkt.discount_factor(2), mkt.discount_factor(4)];
    switch upper(type)
        case 'COC'
            rng(mc.seed);
            V = CoC_pricing_MC(params_AB, scale_factor, mc.N_sim, mc.M, mc.dz, ...
                    mc.N_grid, params_hedge.forward, params_hedge.K1, params_hedge.K2, discount_factors, model, false);

        case 'POP'
            rng(mc.seed);
            V = PoP_pricing_MC(params_AB, scale_factor, mc.N_sim, mc.M, mc.dz, ...
                    mc.N_grid, params_hedge.forward, params_hedge.K1, params_hedge.K2, discount_factors, model, false);

        case 'CHOOSER'
            rng(mc.seed);
            V = Chooser_pricing_MC(params_AB, scale_factor, mc.N_sim, mc.M, mc.dz, ...
                    mc.N_grid, params_hedge.forward, params_hedge.K2, discount_factors, model, false);

        otherwise
            error('price_exotic_AB:badType', ...
                'Unknown exotic "%s". Use ''CoC'', ''PoP'' or ''Chooser''.', type);
    end
end
