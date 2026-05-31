function chooser_price = Chooser_pricing_analytic(params, scale_factor, df, F_t0_T2, K2, diagnostic)
% CHOOSER_PRICING_ANALYTIC Computes the analytic price of a Chooser Option
% under the Minimal Additive (MA) model.
%
% Inputs:
%   params       - 1x2 vector containing the tail decay parameters: [alpha, beta]
%                  where alpha > 0 (left tail) and beta > 0 (right tail).
%   scale_factor - 1x2 vector of scaled volatilities: [sigma_t1, sigma_T2]
%                  where sigma_t1 is at decision time, sigma_T2 is at maturity.
%   df           - 1x2 vector of discount factors from t0: [B(t0, t1), B(t0, T2)]
%   F_t0_T2      - Forward price of the underlying at t0 for maturity T2.
%   K2           - Strike price of the Chooser option at T2 (scalar or vector).
%   diagnostic   - (Optional) Boolean flag for debug prints. Default is false.
%
% Output:
%   chooser_price - The fair value of the Chooser option at t0.

    % =========================================================================
    % DEFAULT ARGUMENT HANDLING
    % =========================================================================
    if nargin < 6 || isempty(diagnostic)
        diagnostic = false;
    end
    
    if diagnostic, tic; end % Start execution timer
    
    % =========================================================================
    % DIAGNOSTICS: INPUT VALIDATION
    % =========================================================================
    if diagnostic
        fprintf('\n==================================================\n');
        fprintf('  RUNNING DIAGNOSTICS FOR Chooser_pricing_analytic\n');
        fprintf('==================================================\n');
        fprintf('INPUT PARAMETERS:\n');
        fprintf('  alpha: %10.4f | beta: %10.4f\n', params(1), params(2));
        fprintf('  sigma_t1:             %10.4f\n', scale_factor(1));
        fprintf('  sigma_T2:             %10.4f\n', scale_factor(2));
        fprintf('  K2 (Chooser Strike):  %10.4f\n', K2(1)); % Indexing prevents vector print overflow
        fprintf('  F(t0,T2) Initial Fwd: %10.4f\n', F_t0_T2);
    end
    
    % =========================================================================
    % CORE PRICING LOGIC
    % =========================================================================
    
    % 1. Calculate the relative strike 'k' (moneyness for additive processes)
    k = K2 - F_t0_T2;
    
    % 2. Calculate the base Vanilla Call price analytically at maturity T2
    % Note: We pass the relative strike 'k', NOT the absolute strike 'K2', 
    % because the analytic function evaluates the zero-mean process.
    call_price_t2 = call_pricing_analytic_MA(k, params, scale_factor(2), df(2), diagnostic);
    
    % 3. Apply Put-Call Parity to derive the base Vanilla Put price at t0
    % P = C + B(t0,T2) * (K - F)
    put_price_t2 = call_price_t2 + df(2) * k;
    
    % 4. Evaluate the Option to Choose (Call on the process at t1)
    % The payoff is evaluated at t1, but we discount directly to t0 using df(2).
    % Why df(2)? Because the choice payoff is multiplied by B(t1, T2) at time t1.
    % Discounting that to t0 gives: B(t0, t1) * B(t1, T2) = B(t0, T2) = df(2).
    price_choice = call_pricing_analytic_MA(k, params, scale_factor(1), df(2), diagnostic);
    
    % 5. Final Chooser Price: Base Put + Option to Choose
    chooser_price = put_price_t2 + price_choice;
    
    % =========================================================================
    % DIAGNOSTICS: OUTPUT REPORTING
    % =========================================================================
    if diagnostic
        elapsed = toc;
        fprintf('\nSTEP 3 - FINAL OUTPUTS:\n');
        fprintf('  Execution Time:    %10.4f seconds\n', elapsed);
        fprintf('  Base Call (T2):    %10.4f\n', call_price_t2(1));
        fprintf('  Base Put (T2):     %10.4f\n', put_price_t2(1));
        fprintf('  Expected Choice:   %10.4f\n', price_choice(1));
        fprintf('  Final Price:       %10.4f\n', chooser_price(1));
        fprintf('==================================================\n\n');
    end
end