function sse = objective_function_MA(x, df, yf, sigma_ATM, mon_mod, c_mkt)
% OBJECTIVE_FUNCTION_MA  SSE between MA model and market call prices.
%   Prices the (M x N) surface via price_MA, sums squared residuals over quoted
%   points only (NaNs ignored), and penalizes model failures.
%
% INPUTS:
%   x         : [alpha, beta]
%   df, yf, sigma_ATM : (M x 1) discount factors, year fractions, ATM vols
%   mon_mod   : (M x N) modified moneyness chi (NaNs allowed)
%   c_mkt     : (M x N) market call prices (NaNs allowed)
% OUTPUT:
%   sse       : scalar SSE to minimize
    alpha = x(1);
    beta  = x(2);
    c_mod = price_MA(x, df, yf, sigma_ATM, mon_mod);

    % Penalty: model returned NaN on a valid market target.
    if any(isnan(c_mod(:)) & ~isnan(c_mkt(:)))
        sse = 1e10 + 1e5 * ((alpha - 1)^2 + (beta - 1)^2);
        return;
    end

    resid_mat = c_mkt - c_mod;
    valid_idx = ~isnan(resid_mat);

    % Penalty: surface collapsed to zero (guards a spurious SSE = 0).
    if sum(valid_idx(:)) == 0 || all(c_mod(valid_idx) == 0)
        sse = 1e10 + 1e5 * ((alpha - 1)^2 + (beta - 1)^2);
        return;
    end

    sse = sum(resid_mat(valid_idx).^2);

end