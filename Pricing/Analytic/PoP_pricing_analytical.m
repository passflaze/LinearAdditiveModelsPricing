function price = PoP_pricing_analytical(params, scale_factor, F_t0_T2, K1, K2, df, diagnostic)
% POP_PRICING_ANALYTICAL Computes the Semi-Analytic price of a Put-on-Put 
% by numerically integrating the exact Put price over the MA density.
%
% Inputs:
%   params            - Tail decay parameters [alpha, beta]
%   scale_factor      - Vector of integrated volatilities [sigma*sqrt(t1), sigma*sqrt(t2)]
%   F_t0_T2           - Initial Forward price
%   K1                - Strike of the Compound Option
%   K2                - Strike of the Underlying Put
%   df                - Vector of discount factors [df_t1, df_t2]
%   diagnostic        - (Optional) Boolean flag to enable debug prints. Default is false.

    % =========================================================================
    % DEFAULT ARGUMENT HANDLING
    % =========================================================================
    if nargin < 7 || isempty(diagnostic)
        diagnostic = false;
    end
    alpha = params(1);
    beta = params(2);
    
    if diagnostic, tic; end % Start execution timer
    
    % =========================================================================
    % DIAGNOSTICS: INPUT VALIDATION
    % =========================================================================
    if diagnostic
        fprintf('\n==================================================\n');
        fprintf('  RUNNING DIAGNOSTICS FOR PoP_pricing_analytical\n');
        fprintf('==================================================\n');
        fprintf('INPUT PARAMETERS:\n');
        fprintf('  alpha: %10.4f | beta: %10.4f\n', alpha, beta);
        fprintf('  K1 (Compound Strike): %10.4f\n', K1);
        fprintf('  K2 (Vanilla Strike):  %10.4f\n', K2);
        fprintf('  F(t0,T2) Initial Fwd: %10.4f\n', F_t0_T2);
        fprintf('  scale_factor(1): %10.4f | df(1): %10.4f | df(2): %10.4f\n', scale_factor(1), df(1), df(2));
        fprintf('\n--- PoP ANALYTICAL ENGINE START ---\n');
    end
    
    % 1. Forward conditional discount factor
    df_t1_t2 = df(2) / df(1);
    
    % 2. Define the inner Put price as a function of the stochastic increment 'x'
    % Utilizing Put-Call Parity: P = C + DF * (K - F)
    put_call_parity_term = @(x) df_t1_t2 * (K2 - F_t0_T2 - x);
    put_val = @(x) call_pricing_analytic_increments_MA(K2 - F_t0_T2 - x, params, scale_factor, df_t1_t2) + put_call_parity_term(x);
    
    % 3. Find the critical increment (x_star) where the Put is exactly ATM (Value = K1)
    objfun = @(x) K1 - put_val(x);
    
    % Capture the exit flag and the residual error (fval) to verify convergence
    [x_star, fval, exitflag] = fzero(objfun, 0);
    
    if diagnostic
        fprintf('ROOT FINDING (fzero):\n');
        fprintf('  x_star found:       %10.4f\n', x_star);
        fprintf('  Put(x_star) value:  %10.4f (Target K1: %.4f)\n', put_val(x_star), K1);
        fprintf('  Residual error:     %10.2e\n', fval);
    end
    
    if exitflag ~= 1 && diagnostic
        warning('Diagnostic Warning: fzero did not converge normally! Exit flag: %d', exitflag);
    end
    
    % 4. Define the Integrand: Payoff(x) * PDF(x)
    integrand = @(x) (K1 - put_val(x)) .* pdf_MA(params, scale_factor(1), x);
    
    % 5. Numerically integrate from the critical point x_star to Infinity 
    % (Put value decreases as x increases, so we exercise when x > x_star)
    integral_val = quadgk(integrand, x_star, Inf);
    
    % 6. Final discounting
    price = df(1) * integral_val;
    
    if diagnostic
        elapsed = toc;
        fprintf('\nINTEGRATION (quadgk):\n');
        fprintf('  Raw Integral (Expected Payoff): %10.4f\n', integral_val);
        fprintf('  Final Discounted Price:         %10.4f\n', price);
        fprintf('  Execution Time:                 %10.4f seconds\n', elapsed);
        fprintf('==================================================\n\n');
    end
    
end