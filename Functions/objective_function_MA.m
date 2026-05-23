function sse = objective_function_MA(x, df, yf, sigma_ATM, mon_mod, c_mkt)
% OBJECTIVE_FUNCTION_MA Computes scalar SSE between model and market prices.
% NaNs in c_mkt are treated as zero contribution to the sum (ignored).

    % Generate the candidate option pricing surface matrix
    c_mod = price_MA(x, df, yf, sigma_ATM, mon_mod);

    % Compute point-by-point errors
    residuals = c_mkt - c_mod;

    % Treat NaNs as zero contribution (ignore them in the sum)
    residuals(isnan(residuals)) = 0;

    % Return scalar SSE
    sse = sum(residuals(:).^2);

end