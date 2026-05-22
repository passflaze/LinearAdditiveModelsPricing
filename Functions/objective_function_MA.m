function sse = objective_function_MA(x, df, yf, sigma_ATM, mon_mod, c_mkt)
% COMPUTE_RESIDUALS Generates the vector of differences between model and market.
% NaNs are computationally filtered out to maintain solver analytical stability.

    if x(1) <= 1e-5 || x(2) <= 1e-5
        sse = 1e10; 
        return;
    end
    % Generate the candidate option pricing surface matrix
    c_mod = price_MA(x, df, yf, sigma_ATM, mon_mod);
    
    % Compute point-by-point errors
    diff = c_mkt - c_mod;
    
    % Filter out empty slots (In-the-Money options filtered previously)
    valid_idx = ~isnan(diff);
    
    % Extract valid data points into a column vector for lsqnonlin
    residuals = diff(valid_idx);

    sse = sum(residuals.^2);
end