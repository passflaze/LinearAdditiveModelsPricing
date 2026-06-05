function [price, CI, ft1, call_price_t1, sigma] = Chooser_pricing_MC(params, scale_factor, N_sim, M, dz, forward, K2, discount_factors, model, opts)
% CHOOSER_PRICING_MC  Computes the Chooser-option price via semi-analytic
%                     Monte Carlo under the MA, GL, or AB model.
%
%   At t1, the increment ft1 = f_{T1,T1} is simulated via CDF inversion (FFT).
%   The forward at t1 is rebuilt with the Lemma-2 rescaling (Forward.pdf),
%       F(T1,T2) = forward + fwd_factor * f_{T1,T1},  fwd_factor = B(0,T1)/B(0,T2).
%   The holder picks max(C_t1, P_t1) of the two inner-option values at t1,
%   which is then discounted to t0. Equivalently this equals
%       B(t0,T2) * E0[ max( E1[(S-K2)+], E1[(K2-S)+] ) ].
%
% INPUTS:
%   params          - (vector) model parameters passed to cf and pricing functions
%   scale_factor    - [scale_t1, scale_t2] scaling factors at t1 and t2
%   N_sim           - (scalar) number of Monte Carlo simulated paths
%   M               - (scalar) grid size exponent for the INNER-CALL Lewis-FFT, N = 2^M
%   dz              - (scalar) step size in the z-domain
%   forward         - (scalar) forward price F(t0, T2)
%   K2              - (scalar) strike of the inner options (paid at t2)
%   discount_factors- [B(t0,t1), B(t0,t2)] discount factors
%   model           - (string) model identifier: 'MA', 'GL', or 'AB'
%   opts            - (optional) struct with fields .verbose and .plot (default false)
%
% OUTPUTS:
%   price           - (scalar) estimated fair value of the Chooser at t0
%   CI              - (1x2) 95% confidence interval [lower, upper]
%   ft1             - (vector) simulated increments at t1
%   call_price_t1   - (vector) inner call prices at t1 for each path
%   sigma           - (scalar) sample standard deviation of the discounted payoff

    % --- Options Initialization ---
    if nargin < 11 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'verbose')
        opts.verbose = false;
    end
    if ~isfield(opts, 'plot')
        opts.plot = false;
    end
    % ------------------------------

    
    % if N_sim == 0
    %     % Accuracy sizing decoupled into size_Nsim_MC and cached on disk: the
    %     % pilot (which also rebuilds the FFT CDF) runs ONCE per (option, model,
    %     % inputs) and is reused. Target = 10 bps market-maker spread (Ex.4).
    %     sig_inputs = [params(:); scale_factor(:); forward; K2; ...
    %                   discount_factors(:); M; dz; N_grid];
    %     N_sim = size_Nsim_MC( ...
    %         @(Np) nth_out(5, @Chooser_pricing_MC, params, scale_factor, Np, M, dz, ...
    %                       N_grid, forward, K2, discount_factors, model, false), ...
    %         sprintf('Chooser_%s', model), sig_inputs, struct('ref', 1));
    % end
    if N_sim == 0
        % Run a pilot simulation to estimate standard deviation (e.g., 1000 paths)
        [~, ~, ~, ~, sigma_est] = Chooser_pricing_MC(params, scale_factor, 1000, M, ...
            dz, forward, K2, discount_factors, model);
        target_error = 10 * 1e-4; 
        N_sim = min(ceil(((1.96 * sigma_est) / target_error)^2), 5e7);
        
        if opts.verbose
            fprintf('--- PILOT SIMULATION ---\n');
            fprintf('Estimated Std Dev: %.4f\n', sigma_est);
            fprintf('Target Error:      10 bps (%.4f)\n', target_error);
            fprintf('Required N_sim:    %d\n', N_sim);
            fprintf('------------------------\n');
        end
    end
    

    % Lemma 2 (Forward.pdf): F(T1,T2) = forward + fwd_factor * f_{T1,T1}.
    fwd_factor = discount_factors(1) / discount_factors(2);

    if opts.verbose
        fprintf('\n=================================================================\n');
        fprintf('  CHOOSER PRICING MC — Model: %s\n', model);
        fprintf('=================================================================\n');
        fprintf('  Parameters:       [%.6f, %.6f]\n', params(1), params(2));
        fprintf('  Scale factors:    scale_t1 = %.6f  |  scale_t2 = %.6f\n', scale_factor(1), scale_factor(2));
        fprintf('  Discount factors: B(t0,t1) = %.6f  |  B(t0,t2) = %.6f\n', discount_factors(1), discount_factors(2));
        fprintf('  fwd_factor (L2):  %.6f\n', fwd_factor);
        fprintf('  Forward F(t0,t2): %.6f\n', forward);
        fprintf('  Strike K2: %.6f\n', K2);
        fprintf('  N_sim: %d  |  M: %d  |  dz: %.4f\n', N_sim, M, dz);
        fprintf('-----------------------------------------------------------------\n');
        fprintf('  [1/3] Simulating ft1...\n');
    end

    % --- Simulate the increment ft1 = f_{T1,T1} via CDF inversion ---
    switch model
        case 'MA'
            % params(1) = alpha, params(2) = beta
            ps_plus = params(2)/scale_factor(1); 
            ps_minus = params(1)/scale_factor(1);
            gamma_MA = (1 / params(1)) - (1 /params(2));
            drift_0_t1 = gamma_MA * scale_factor(1);
            
            if opts.verbose
                fprintf('        drift_0_t1 = %.6f  |  std_T1 = %.6f\n', drift_0_t1, std_T1);
            end
            
            ft1 = FA_simulation(N_sim, M, dz, drift_0_t1, ...
                             ps_plus, ps_minus, 0, 0, 1, 'infinite', 1, params, scale_factor(1), opts.plot);
        case 'GL'
            % params(1) = alpha, params(2) = beta
            [cdf_fT1, z_grid]    = lewis_FFT_digital(@cf_GL, M, dz, params, ...
                             scale_factor(1), 1, 'GL', 1, opts.plot);
                             
            if opts.verbose
                fprintf('        z_grid: [%.4f, %.4f]  (%d points)\n', z_grid(1), z_grid(end), numel(z_grid));
            end
            
            ft1        = simulate_from_cdf(cdf_fT1, z_grid, 1, N_sim);
        case 'AB'
            % params(1) = k (kappa), params(2) = eta
            
            [cdf_fT1, z_grid]    = lewis_FFT_digital(@cf_AB, M, dz, params, ...
                             scale_factor(1), 1, 'AB', 1, opts.plot);
                             
            if opts.verbose
                fprintf('        z_grid: [%.4f, %.4f]  (%d points)\n', z_grid(1), z_grid(end), numel(z_grid));
            end
            
            ft1        = simulate_from_cdf(cdf_fT1, z_grid, 1, N_sim);
    end
    
    if opts.verbose
        fprintf('        ft1 — mean: %.6f  |  std: %.6f  |  min: %.4f  |  max: %.4f\n', ...
            mean(ft1), std(ft1), min(ft1), max(ft1));
    end

    % --- Forward price at t1 for each simulated path (Lemma 2 rescaling) ---
    F_t1_T2 = forward + fwd_factor * ft1;
    
    % --- Analytically price the inner options at t1 for each path ---
    if opts.verbose
        fprintf('-----------------------------------------------------------------\n');
        fprintf('  [2/3] Pricing inner options at t1...\n');
    end
    
    % strikes are the dollar-moneyness K2 - F_t1_T2
    strikes      = K2 - F_t1_T2;
    df_t1_t2     = discount_factors(2) / discount_factors(1);
    
    if opts.verbose
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
    
    % Inner put value at t1 via put-call parity: P = C + B(t1,t2)*(K2 - F)
    put_price_t1   = call_price_t1 + df_t1_t2 * strikes;
    
    if opts.verbose
        fprintf('        call_t1 — mean: %.6f  |  std: %.6f\n', mean(call_price_t1), std(call_price_t1));
        fprintf('        put_t1  — mean: %.6f  |  std: %.6f\n', mean(put_price_t1), std(put_price_t1));
        fprintf('-----------------------------------------------------------------\n');
        fprintf('  [3/3] Computing compound payoff and discounting to t0...\n');
    end
    
    payoffs = max(call_price_t1, put_price_t1);
    discounted_payoffs = discount_factors(1) * payoffs;
    
  
    [price, sigma, CI, ~] = normfit(discounted_payoffs);
    
    if opts.verbose
        fprintf('=================================================================\n');
        fprintf('  RESULTS — Model: %s\n', model);
        fprintf('  Price:  %.6f\n', price);
        fprintf('  95%% CI: [%.6f, %.6f]\n', CI(1), CI(2));
        fprintf('  std:  %.6f\n', sigma);
        fprintf('=================================================================\n\n');
    end
end