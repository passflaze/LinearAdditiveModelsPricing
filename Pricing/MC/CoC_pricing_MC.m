function [price, CI, ft1, call_price_t1] = CoC_pricing_MC(params, scale_factor, N_sim, M, dz, N_grid, forward, K1, K2, discount_factors, model, diagnostics)
% COC_PRICING_MC  Computes the Call-on-Call price via semi-analytic
%                 Monte Carlo under the MA, GL, or AB model.
%
%   At t1, the increment ft1 = f_{T1,T1} is simulated via CDF inversion (FFT).
%   The forward at t1 is rebuilt with the Lemma-2 rescaling (Forward.pdf),
%       F(T1,T2) = forward + fwd_factor * f_{T1,T1},  fwd_factor = B(0,T1)/B(0,T2),
%   the inner call value at t1 is then computed analytically for each path and
%   the compound payoff max(C_t1 - K1, 0) is discounted to t0.
%
% INPUTS:
%   params          - (vector) model parameters passed to cf and pricing functions
%   scale_factor    - [scale_t1, scale_t2] scaling factors at t1 and t2
%   N_sim           - (scalar) number of Monte Carlo simulated paths
%   M               - (scalar) grid size exponent, N = 2^M
%   dz              - (scalar) step size in the z-domain
%   N_grid          - (scalar) number of points for the MA CDF grid
%   forward         - (scalar) forward price F(t0, T2)
%   K1              - (scalar) strike of the compound option (paid at t1)
%   K2              - (scalar) strike of the inner vanilla call (paid at t2)
%   discount_factors- [B(t0,t1), B(t0,t2)] discount factors
%   model           - (string) model identifier: 'MA', 'GL', or 'AB'
%   diagnostics     - (optional, default true) print the per-model diagnostics.
%
% OUTPUTS:
%   price           - (scalar) estimated fair value of the Call-on-Call at t0
%   CI              - (1x2) 95% confidence interval [lower, upper]
%   ft1             - (vector) simulated increments at t1
%   call_price_t1   - (vector) inner call prices at t1 for each path

    if nargin < 12 || isempty(diagnostics), diagnostics = false; end
    
    % Lemma 2 (Forward.pdf): F(T1,T2) = forward + fwd_factor * f_{T1,T1}.
    fwd_factor = discount_factors(1) / discount_factors(2);
    
    if diagnostics
        fprintf('\n=================================================================\n');
        fprintf('  COC PRICING MC — Model: %s\n', model);
        fprintf('=================================================================\n');
        fprintf('  Parameters:       [%.6f, %.6f]\n', params(1), params(2));
        fprintf('  Scale factors:    scale_t1 = %.6f  |  scale_t2 = %.6f\n', scale_factor(1), scale_factor(2));
        fprintf('  Discount factors: B(t0,t1) = %.6f  |  B(t0,t2) = %.6f\n', discount_factors(1), discount_factors(2));
        fprintf('  fwd_factor (L2):  %.6f\n', fwd_factor);
        fprintf('  Forward F(t0,t2): %.6f\n', forward);
        fprintf('  Strike K1: %.6f  |  Strike K2: %.6f\n', K1, K2);
        fprintf('  N_sim: %d  |  M: %d  |  dz: %.4f\n', N_sim, M, dz);
        fprintf('-----------------------------------------------------------------\n');
    end
    
    % --- Simulate the increment ft1 = f_{T1,T1} via CDF inversion ---
    if diagnostics, fprintf('  [1/3] Simulating ft1...\n'); end
    
    switch model
        case 'MA'
            % params(1) = alpha, params(2) = beta
            ps_plus    = params(2) / scale_factor(1);
            ps_minus   = params(1) / scale_factor(1);
            std_T1     = (1/params(1)^2) + (1/params(2)^2);
            z_grid_std = linspace(-100 * std_T1, 100 * std_T1, N_grid)';
            gamma_MA   = (1/params(1)) - (1/params(2));
            drift_0_t1 = gamma_MA * scale_factor(1);
            
            if diagnostics
                fprintf('        drift_0_t1 = %.6f  |  std_T1 = %.6f\n', drift_0_t1, std_T1);
                fprintf('        z_grid_std: [%.4f, %.4f]  (%d points)\n', z_grid_std(1), z_grid_std(end), N_grid);
            end
            
            ft1 = FA_simulation(N_sim, M, dz, drift_0_t1, ...
                             ps_plus, ps_minus, 0, 0, 1, 'infinite', 1, params, scale_factor(1), z_grid_std);
        case 'GL'
            % params(1) = alpha, params(2) = beta
            [cdf_fT1, z_grid] = lewis_FFT_digital(@cf_GL, M, dz, params, ...
                                    scale_factor(1), 1, 'GL', 1);
            
            if diagnostics, fprintf('        z_grid: [%.4f, %.4f]  (%d points)\n', z_grid(1), z_grid(end), numel(z_grid)); end
            
            ft1 = simulate_from_cdf(cdf_fT1, z_grid, 1, N_sim);
        case 'AB'
            % params(1) = k (kappa), params(2) = eta
            [cdf_fT1, z_grid] = lewis_FFT_digital(@cf_AB, M, dz, params, ...
                                    scale_factor(1), 1, 'AB', 1);
            
            if diagnostics, fprintf('        z_grid: [%.4f, %.4f]  (%d points)\n', z_grid(1), z_grid(end), numel(z_grid)); end
            
            ft1 = simulate_from_cdf(cdf_fT1, z_grid, 1, N_sim);
    end
    
    if diagnostics
        fprintf('        ft1 — mean: %.6f  |  std: %.6f  |  min: %.4f  |  max: %.4f\n', ...
            mean(ft1), std(ft1), min(ft1), max(ft1));
    end
    
    % --- Forward price at t1 for each simulated path (Lemma 2 rescaling) ---
    F_t1_T2 = forward + fwd_factor * ft1;
    
    % --- Analytically price the inner call at t1 for each path ---
    if diagnostics
        fprintf('-----------------------------------------------------------------\n');
        fprintf('  [2/3] Pricing inner call at t1...\n');
    end
    
    strikes  = K2 - F_t1_T2;
    df_t1_t2 = discount_factors(2) / discount_factors(1);
    
    if diagnostics
        fprintf('        df(t1,t2) = %.6f\n', df_t1_t2);
        fprintf('        strikes — mean: %.6f  |  std: %.6f\n', mean(strikes), std(strikes));
    end
    
    switch model
        case 'MA'
            call_price_t1 = df_t1_t2 * lewis_FFT_call(@cf_MA_FA, M, dz, params, ...
                                scale_factor, strikes, 1, 'MA', fwd_factor);
        case 'GL'
            call_price_t1 = df_t1_t2 * lewis_FFT_call(@cf_increment_GL, M, dz, params, ...
                                scale_factor, strikes, 1, 'GL', fwd_factor);
        case 'AB'
            call_price_t1 = df_t1_t2 * lewis_FFT_call(@cf_increment_AB, M, dz, params, ...
                                scale_factor, strikes, 1, 'AB', fwd_factor);
    end
    
    if diagnostics
        fprintf('        call_t1 — mean: %.6f  |  std: %.6f\n', mean(call_price_t1), std(call_price_t1));
    end
    
    % --- Compound payoff at t1: max(C_t1 - K1, 0) ---
    if diagnostics
        fprintf('-----------------------------------------------------------------\n');
        fprintf('  [3/3] Computing compound payoff and discounting to t0...\n');
    end
    
    compound_payoff_t1 = max(call_price_t1 - K1, 0);
    discounted_payoffs = discount_factors(1) * compound_payoff_t1;
    
    price   = mean(discounted_payoffs);
    std_err = std(discounted_payoffs) / sqrt(N_sim);
    CI      = [price - 1.96 * std_err, price + 1.96 * std_err];
    
    if diagnostics
        fprintf('=================================================================\n');
        fprintf('  RESULTS — Model: %s\n', model);
        fprintf('  Price:  %.6f\n', price);
        fprintf('  95%% CI: [%.6f, %.6f]\n', CI(1), CI(2));
        fprintf('=================================================================\n\n');
    end
end