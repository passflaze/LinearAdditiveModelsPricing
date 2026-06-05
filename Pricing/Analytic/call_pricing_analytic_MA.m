function call_price = call_pricing_analytic_MA(strike,params, scale_factor, df, diagnostic)
% CALL_PRICING_ANALYTIC_MA Computes the exact analytic price of a Vanilla Call
% evaluated under the Infinite Activity (Marginal f_t) Minimal Additive model.
%
% INPUTS:
%   strike       - The effective strike (can be a vector)
%   params       - [alpha, beta] left/right tail decay parameters
%   scale_factor - Integrated volatility at time t: sigma*sqrt(t) (scalar)
%   df           - Discount factor for time t
%   diagnostic   - (Optional) Boolean flag for debug prints. Default is false.
% OUTPUT:
%   call_price   - The discounted analytic price of the Call

    if nargin < 6 || isempty(diagnostic)
        diagnostic = false;
    end
    alpha = params(1);
    beta=params(2);
    scale_factor = max(scale_factor);
    % 1. Shift calculation (mu_t)
    gamma_MA = (1 / alpha) - (1 / beta);
    mu_t = gamma_MA * scale_factor;
    
    % --- DIAGNOSTIC BLOCK ---
    if diagnostic
        fprintf('  [Marginal Call Analytic] Eval on %d points | mu_t: %.4f | strike range: [%.4f, %.4f]\n', ...
                length(strike), mu_t, min(strike), max(strike));
    end
    
    % 2. Tail parameters for the time t
    pt_plus  = beta  / scale_factor;
    pt_minus = alpha / scale_factor;
    
    % 3. Scaling coefficients based on the Infinite Activity formula
    denominator = pt_plus + pt_minus;
    coeff_left  = pt_plus  / (pt_minus * denominator);
    coeff_right = pt_minus / (pt_plus  * denominator);
    
    % 4. Initialize the output array to match the size of the input strikes
    call_price_undiscounted = zeros(size(strike));
    
    % 5. Create logical masks for the piecewise integration
    idx_right = strike >= mu_t;
    idx_left  = strike <  mu_t; 
    
    % 6. Evaluate the analytic integrals based on the Strike Regime
    
    % Regime 1: k >= mu_t (Right side)
    if any(idx_right)
        call_price_undiscounted(idx_right) = coeff_right .* exp(-pt_plus .* (strike(idx_right) - mu_t));
    end
    
    % Regime 2: k < mu_t (Left side)
    if any(idx_left)
        call_price_undiscounted(idx_left) = -strike(idx_left) + coeff_left .* exp(pt_minus .* (strike(idx_left) - mu_t));
    end
    
    % 7. Apply the discount factor
    call_price = df * call_price_undiscounted;
    
end