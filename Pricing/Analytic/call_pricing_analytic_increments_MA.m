function call_price_t1 = call_pricing_analytic_increments_MA(strike,params, scale_factor, df, diagnostic)
% CALL_PRICING_ANALYTIC Computes the exact analytic price of a Vanilla Call 
% evaluated at t1 (expiring at t2) under the Minimal Additive model.
%
% Inputs:
%   strike - The effective strike evaluated at t1 (can be a vector)
%   alpha      - Left tail decay parameter
%   beta       - Right tail decay parameter
%   scale_factor     - Vector of integrated volatilities [sigma*sqrt(t1), sigma*sqrt(t2)]
%   df   - Discount factor from t1 to t2
%   diagnostic - (Optional) Boolean flag for debug prints. Default is false.
%
% Output:
%   call_price_t1 - The discounted analytic price of the Call in t1

    % =========================================================================
    % DEFAULT ARGUMENT HANDLING
    % =========================================================================
    if nargin < 5 || isempty(diagnostic)
        diagnostic = false;
    end
    alpha = params(1);
    beta = params(2);
    % 1. Shift calculation (deltamu)
    gamma_MA = (1 / alpha) - (1 / beta);
    deltamu = gamma_MA * (scale_factor(2) - scale_factor(1));
    
    % --- DIAGNOSTIC BLOCK ---
    % Kept minimal to avoid terminal flooding during quadgk/fzero iterations
    if diagnostic
        fprintf('  [Inner Call Analytic] Eval on %d points | deltamu: %.4f | strike range: [%.4f, %.4f]\n', ...
                length(strike), deltamu, min(strike), max(strike));
    end
    
    % 2. Tail parameters for the total time (t) and intermediate time (s)
    pt_plus  = beta  / scale_factor(2);
    ps_plus  = beta  / scale_factor(1);
    pt_minus = alpha / scale_factor(2);
    ps_minus = alpha / scale_factor(1);
    
    % Point mass size 'c' (probability of zero jumps)
    c = (pt_plus * pt_minus) / (ps_plus * ps_minus);
    
    % Scaling coefficients for the continuous density branches
    denominator = pt_plus + pt_minus;
    A_plus = c * ((ps_plus - pt_plus) * (ps_minus + pt_plus)) / denominator;
    A_minus  = c * ((ps_minus - pt_minus) * (ps_plus + pt_minus)) / denominator;
    
    % 4. Initialize the output array to match the size of the input strikes
    call_price = zeros(size(strike));
    
    % 5. Create logical masks for the piecewise integration
    idx_right = strike >= deltamu;
    idx_left  = strike <  deltamu; 
    
    % 6. Evaluate the analytic integrals
    % Regime 1: Strike is to the right of the jump (Uses pt_plus)
    if any(idx_right)
        call_price(idx_right) = (A_plus / pt_plus^2) .* exp(-pt_plus .* (strike(idx_right) - deltamu));
    end
    
    % Regime 2: Strike is to the left of the jump (Uses pt_minus and the linear -k term)
    if any(idx_left)
        call_price(idx_left) = -strike(idx_left) + (A_minus / pt_minus^2) .* exp(-pt_minus .* (deltamu - strike(idx_left)));
    end
    
    % 7. Apply the discount factor for the t1 -> t2 period
    call_price_t1 = df * call_price;
    
end