function price = CoC_pricing_FULLanalytical(params, scale_factor, F_t0_T2, K1, K2, df, diagnostic)
% COC_PRICING_FULLANALYTICAL Computes the exact analytic price of a Call-on-Call
% by applying the closed-form primitives of the MA density (Section 5.3.2).
%
% INPUTS:
%   params       - Tail decay parameters [alpha, beta]
%   scale_factor - Vector of integrated volatilities [sigma*sqrt(t1), sigma*sqrt(t2)]
%   F_t0_T2      - Initial forward price
%   K1           - Strike of the compound option
%   K2           - Strike of the underlying call
%   df           - Vector of discount factors [df_t1, df_t2]
%   diagnostic   - (Optional) Boolean flag to enable debug prints. Default is false.
% OUTPUT:
%   price        - The discounted analytic price of the Call-on-Call

    if nargin < 7 || isempty(diagnostic)
        diagnostic = false;
    end
    alpha = params(1);
    beta = params(2);
    
    if diagnostic, tic; end

    % 1. Forward conditional discount factor
    df_t1_t2 = df(2) / df(1);
    
    % 2. Define the inner Call price as a function of the stochastic increment 'x'
    call_val = @(x) call_pricing_analytic_increments_MA(K2 - F_t0_T2 - x, params, scale_factor, df_t1_t2);
    
    % 3. Find the critical increment (x_star) where the Call is exactly ATM (Value = K1)
    objfun = @(x) call_val(x) - K1;
    [x_star, fval, exitflag] = fzero(objfun, 0);
    
    if diagnostic
        fprintf('ROOT FINDING (fzero):\n');
        fprintf('  x_star found:       %10.4f\n', x_star);
        fprintf('  Call(x_star) value: %10.4f (Target K1: %.4f)\n', call_val(x_star), K1);
        fprintf('  Residual error:     %10.2e\n', fval);
    end
    
    if exitflag ~= 1 && diagnostic
        warning('Diagnostic Warning: fzero did not converge normally! Exit flag: %d', exitflag);
    end
    
    % 4. Fully analytical integration (Section 5.3.2)

    % --- 4.1 Define Process and Density Parameters ---
    s1 = scale_factor(1);
    s2 = scale_factor(2);
    
    p1_minus = alpha / s1;
    p1_plus  = beta / s1;
    p2_minus = alpha / s2;
    p2_plus  = beta / s2;
    
    gamma_val = (1/alpha) - (1/beta);
    mu_T1 = gamma_val * s1;
    mu_T2 = gamma_val * s2;
    delta_mu = mu_T2 - mu_T1;
    
    x_hat = K2 - F_t0_T2 - delta_mu;
    
    % --- 4.2 Define Constants for the Inner Call (Eq 44-47) ---
    c_val   = (p2_plus * p2_minus) / (p1_plus * p1_minus);
    A_minus = c_val * ((p1_minus - p2_minus) * (p1_plus + p2_minus)) / (p2_plus + p2_minus);
    A_plus  = c_val * ((p1_plus - p2_plus) * (p1_minus + p2_plus)) / (p2_plus + p2_minus);
    
    C_out = df_t1_t2 * (A_plus / p2_plus^2) * exp(-p2_plus * x_hat);
    C_in  = df_t1_t2 * (A_minus / p2_minus^2) * exp(p2_minus * x_hat);
    C_lin = df_t1_t2 * (K2 - F_t0_T2);
    
    % --- 4.3 Define Constants for the Marginal Density (Eq 48) ---
    const_D = (p1_plus * p1_minus) / (p1_plus + p1_minus);
    D_L = const_D * exp(-p1_minus * mu_T1);
    D_R = const_D * exp(p1_plus * mu_T1);
    
    % --- 4.4 Define the Analytical Primitives (Eq 49-52) ---
    % Regime 1: x <= x_hat AND x < mu_T1
    Psi1 = @(x) D_L * ( (C_out / (p1_minus + p2_plus)) * exp((p1_minus + p2_plus) * x) - (K1 / p1_minus) * exp(p1_minus * x) );
    
    % Regime 2: x <= x_hat AND x >= mu_T1
    Psi2 = @(x) D_R * ( (C_out / (p2_plus - p1_plus)) * exp((p2_plus - p1_plus) * x) + (K1 / p1_plus) * exp(-p1_plus * x) );
    
    % Regime 3: x > x_hat AND x < mu_T1
    Psi3 = @(x) D_L * ( df_t1_t2 * ((p1_minus * x - 1) / p1_minus^2) .* exp(p1_minus * x) ...
                      - ((C_lin + K1) / p1_minus) * exp(p1_minus * x) ...
                      + (C_in / (p1_minus - p2_minus)) * exp((p1_minus - p2_minus) * x) );
                      
    % Regime 4: x > x_hat AND x >= mu_T1
    Psi4 = @(x) D_R * ( -df_t1_t2 * ((p1_plus * x + 1) / p1_plus^2) .* exp(-p1_plus * x) ...
                      + ((C_lin + K1) / p1_plus) * exp(-p1_plus * x) ...
                      - (C_in / (p1_plus + p2_minus)) * exp(-(p1_plus + p2_minus) * x) );
                      
    % --- 4.5 Domain Partitioning ---
    % Identify structural breaks and sort intervals strictly >= x_star
    all_breaks = [x_star, x_hat, mu_T1];
    bounds = sort(all_breaks(all_breaks >= x_star));
    
    % Ensure x_star is always the leading boundary
    if isempty(bounds) || bounds(1) ~= x_star
        bounds = [x_star, bounds];
    end
    
    bounds = unique(bounds); % Remove overlapping boundaries
    bounds = [bounds, Inf];  % Close the domain to infinity
    
    % --- 4.6 Exact Integration Loop ---
    integral_val = 0;
    
    for i = 1:(length(bounds)-1)
        a = bounds(i);
        b = bounds(i+1);
        
        if b == a
            continue;
        end
        
        % Pick a midpoint to identify the active regime for this segment
        if isinf(b)
            mid = a + 1; 
        else
            mid = (a + b) / 2;
        end
        
        % Evaluate the active primitive [Psi(b) - Psi(a)]
        % Note: Psi4(+Inf) is mathematically 0 because exponents are negative.
        if mid <= x_hat && mid < mu_T1
            val_a = Psi1(a);
            if isinf(b), val_b = 0; else, val_b = Psi1(b); end
            
        elseif mid <= x_hat && mid >= mu_T1
            val_a = Psi2(a);
            if isinf(b), val_b = 0; else, val_b = Psi2(b); end
            
        elseif mid > x_hat && mid < mu_T1
            val_a = Psi3(a);
            if isinf(b), val_b = 0; else, val_b = Psi3(b); end
            
        elseif mid > x_hat && mid >= mu_T1
            val_a = Psi4(a);
            if isinf(b), val_b = 0; else, val_b = Psi4(b); end
        end
        
        integral_val = integral_val + (val_b - val_a);
    end
    
    % 5. Final discounting
    price = df(1) * integral_val;
    
    if diagnostic
        elapsed = toc;
        fprintf('\nINTEGRATION (Exact Primitives):\n');
        fprintf('  Raw Integral (Expected Payoff): %10.4f\n', integral_val);
        fprintf('  Final Discounted Price:         %10.4f\n', price);
        fprintf('  Execution Time:                 %10.4f seconds\n', elapsed);
        fprintf('==================================================\n\n');
    end
    
end