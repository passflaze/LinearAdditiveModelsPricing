function g = price_hedge_basket(basket, params_AB, scale_factor_vanilla, mkt, params_hedge, mc, bumps)
% PRICE_HEDGE_BASKET  Price + Greeks of a basket of vanilla/future hedging
%                     instruments, in one call.
%
%   Thin loop over greeks_vanilla_AB so that run_ex6 and hedge_backtest share
%   exactly the same per-instrument pricing/Greek conventions.
%
% INPUTS:
%   basket               - (N x 1 struct array) specs with .kind .K .mat
%   params_AB            - (2x1) [k; eta]
%   scale_factor_vanilla - (vector) sigma_t .* sqrt(yf) per maturity
%   mkt                  - (struct) market data
%   params_hedge         - (struct) hedging parameters (forward, ...)
%   mc, bumps            - (struct) numerical / bump settings
%
% OUTPUT:
%   g - (N x 1 struct array) with .kind .K .mat .price .delta .gamma .vega

    N = numel(basket);
    g = repmat(struct('kind', '', 'K', NaN, 'mat', NaN, ...
                      'price', 0, 'delta', 0, 'gamma', 0, 'vega', 0), N, 1);

    for j = 1:N
        % evalc keeps the console clean (suppresses any calibration penalties)
        g(j) = greeks_vanilla_AB(basket(j).kind, basket(j).K, basket(j).mat, ...
                           params_AB, scale_factor_vanilla, mkt, params_hedge, mc, bumps);
    end
end
