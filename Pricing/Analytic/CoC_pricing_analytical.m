function price = CoC_pricing_analytical(params, scale_factor, F_t0_T2, K1, K2, df, diagnostic)
% COC_PRICING_ANALYTICAL Computes the Semi-Analytic price of a Call-on-Call
% by numerically integrating the exact Call price over the MA density.
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
%   price        - The discounted semi-analytic price of the Call-on-Call

    if nargin < 8 || isempty(diagnostic)
        diagnostic = false;
    end
    alpha = params(1);
    beta = params(2);

    if diagnostic, tic; end

    if diagnostic
        fprintf('\n==================================================\n');
        fprintf('  RUNNING DIAGNOSTICS FOR CoC_pricing_analytical\n');
        fprintf('==================================================\n');
        fprintf('INPUT PARAMETERS:\n');
        fprintf('  alpha: %10.4f | beta: %10.4f\n', alpha, beta);
        fprintf('  K1 (Compound Strike): %10.4f\n', K1);
        fprintf('  K2 (Vanilla Strike):  %10.4f\n', K2);
        fprintf('  F(t0,T2) Initial Fwd: %10.4f\n', F_t0_T2);
        fprintf('  scale_factor(1): %10.4f | df(1): %10.4f | df(2): %10.4f\n', scale_factor(1), df(1), df(2));
        fprintf('\n--- CoC ANALYTICAL ENGINE START ---\n');
    end

    % 1. Forward conditional discount factor and Lemma-2 forward rescaling.
    %    F(T1,T2) = F(t0,T2) + fwd_factor * f_{T1,T1},  fwd_factor = B(0,T1)/B(0,T2).
    df_t1_t2   = df(2) / df(1);
    fwd_factor = df(1) / df(2);

    % 2. Define the inner Call price as a function of the increment x = f_{T1,T1}.
    %    The inner-call increment law carries fwd_factor exactly as in cf_MA_FA;
    %    since cf_MA_FA(u,[s1,s2],fwd) == cf_MA_FA(u,[fwd*s1,s2],1), the closed
    %    form is obtained by rescaling the s-leg scale by fwd_factor. The dollar
    %    moneyness seen at t1 is K2 - F_t1_T2 = K2 - F(t0,T2) - fwd_factor*x.
    scale_inc = [fwd_factor * scale_factor(1), scale_factor(2)];
    call_val  = @(x) call_pricing_analytic_increments_MA(K2 - F_t0_T2 - fwd_factor*x, params, scale_inc, df_t1_t2);
    
    % 3. Find the critical increment (x_star) where the Call is exactly ATM (Value = K1)
    objfun = @(x) call_val(x) - K1;
    
    % Capture the exit flag and the residual error (fval) to verify convergence
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
    
    % 4. Define the Integrand: Payoff(x) * PDF(x)
    integrand = @(x) (call_val(x) - K1) .* pdf_MA(params, scale_factor(1), x);
    
    % 5. Numerically integrate from the critical point x_star to Infinity 
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