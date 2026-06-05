function comparison_function_increments_MA(params, scale_factor, df, M, dz,K2,fwd)
% COMPARISON_FUNCTION_INCREMENTS_MA Compares the analytic and FFT call prices
% for the MA increment (t1 -> t2) to verify their consistency.
%
% INPUTS:
%   params       - [alpha, beta] left/right tail decay parameters
%   scale_factor - Vector of integrated volatilities [sigma*sqrt(t1), sigma*sqrt(t2)]
%   df           - Vector of discount factors [B(t0,t1), B(t0,t2)]
%   M            - Controls FFT grid size (N = 2^M)
%   dz           - Step size for the frequency grid
%   K2           - Strike of the inner call (paid at t2)
%   fwd          - Forward price F(t1, T2)
% OUTPUT:
%   none (prints a comparison table of analytic vs FFT prices)

    % 1. Derived Tail Parameters
    alpha = params(1); beta = params(2);
    pt_plus  = beta  / scale_factor(2);
    pt_minus = alpha / scale_factor(2);
    ps_plus  = beta  / scale_factor(1);
    ps_minus = alpha / scale_factor(1);
    gamma_MA = (1/alpha) - (1/beta);
    drift_0_t1 = gamma_MA*(scale_factor(1));
    N_sim = 1e6;
    df_t1_t2 = df(2)/df(1);
    
    ft1 = FA_simulation(N_sim, M, dz,drift_0_t1,pt_plus, pt_minus ,...
                             ps_plus, ps_minus, 1, 'infinite', 1, params,scale_factor);
    F_t1_T2 = fwd + ft1;
 
    strikes = K2 - F_t1_T2;
    % 2. Analytic Price
    fprintf('Running Analytic Pricing...\n');
    price_analytic = call_pricing_analytic_increments_MA(strikes, params, scale_factor, df_t1_t2, false);

    % 3. FFT Price
    fprintf('Running FFT Pricing...\n');
    % Call the FFT pricing function, applying the discount factor to match analytic
    price_fft_raw = lewis_FFT_call(@cf_MA_FA, M, dz, params, scale_factor, strikes,1, 'MA');
    price_fft = df_t1_t2 * price_fft_raw;

    % 4. Comparison & Diagnostics
    fprintf('\n--- COMPARISON RESULTS ---\n');
    
    % Ensure column vectors for clean table output
    strikes = strikes(:);
    price_analytic = price_analytic(:);
    price_fft = price_fft(:);
    
    diff_abs = abs(price_analytic - price_fft);
    diff_bps = diff_abs * 10000;

    fprintf('%-10s | %-15s | %-15s | %-12s | %-10s\n', 'Strike', 'Analytic', 'FFT', 'Diff (Abs)', 'Diff (bps)');
    fprintf('----------------------------------------------------------------------\n');
    
    % Display up to 20 samples to avoid terminal flooding
    num_display = min(length(strikes), 20);
    
    % Pick indices evenly spaced across the strike vector
    idx_display = round(linspace(1, length(strikes), num_display));
    
    for i = 1:num_display
        idx = idx_display(i);
        fprintf('%10.4f | %15.6f | %15.6f | %12.6f | %10.2f\n', ...
            strikes(idx), price_analytic(idx), price_fft(idx), diff_abs(idx), diff_bps(idx));
    end
    
    fprintf('----------------------------------------------------------------------\n');
    fprintf('Max Absolute Error: %e\n', max(diff_abs));
    fprintf('Mean Absolute Error: %e\n', mean(diff_abs));
    fprintf('----------------------------------------------------------------------\n');

end