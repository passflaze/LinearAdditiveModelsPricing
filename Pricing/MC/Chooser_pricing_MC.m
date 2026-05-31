function [price, CI, ft1, call_price_t1] = Chooser_pricing_MC(params, scale_factor, N_sim, M, dz, N_grid, forward, K2, discount_factors, model)
% COC_PRICING_MC_PROVA  Computes the Call-on-Call price via semi-analytic
%                       Monte Carlo under the MA, GL, or AB model.
%
%   At t1, the increment ft1 is simulated via CDF inversion (FFT). The
%   inner call price at t1 is then computed analytically for each path.
%   The compound payoff max(C_t1 - K1, 0) is discounted to t0.
%
% INPUTS:
%   params          - (vector) model parameters passed to cf and pricing functions
%   scale_factor    - [scale_t1, scale_t2] scaling factors at t1 and t2
%   N_sim           - (scalar) number of Monte Carlo simulated paths
%   M               - (scalar) grid size exponent, N = 2^M
%   dz              - (scalar) step size in the z-domain
%   N_grid          - (scalar) number of points for the CDF grid
%   forward         - (scalar) forward price F(t0, T2)
%   K1              - (scalar) strike of the compound option (paid at t1)
%   K2              - (scalar) strike of the inner vanilla call (paid at t2)
%   discount_factors- [B(t0,t1), B(t0,t2)] discount factors
%   model           - (string) model identifier: 'MA', 'GL', or 'AB'
%
% OUTPUTS:
%   price           - (scalar) estimated fair value of the Call-on-Call at t0
%   CI              - (1x2) 95% confidence interval [lower, upper]
%   ft1             - (vector) simulated increments at t1
%   call_price_t1   - (vector) inner call prices at t1 for each path

    % --- Simulate the increment ft1 via CDF inversion ---
    switch model
        case 'MA'
            % params(1) = alpha, params(2) = beta
            ps_plus = params(2)/scale_factor(1); ps_minus = params(1)/scale_factor(1);
            std_T1     = (1/params(1)^2) + (1/params(2)^2);
            z_grid_std = linspace(-100 * std_T1, 100 * std_T1, N_grid)';
            gamma_MA = (1 / params(1)) - (1 /params(2));
            drift_0_t1 = gamma_MA * (scale_factor(1) -0);
            ft1 = FA_simulation(N_sim, M, dz, drift_0_t1, ...
                             ps_plus, ps_minus, 0, 0, 1, 'infinite',1, params, scale_factor(1), z_grid_std);

        case 'GL'
            % params(1) = alpha, params(2) = beta
            % std_T1     = scale_factor(1) * sqrt(psi(1, params(1)) + psi(1, params(2)));
            % z_grid_std = linspace(-10 * std_T1, 10 * std_T1, N_grid)';
            [cdf_fT1, z_grid]    = lewis_FFT_digital(@cf_GL, M, dz, params, ...
                             scale_factor(1), 1, 'GL', 1);
            ft1        = simulate_from_cdf(cdf_fT1, z_grid, 1, N_sim);

        case 'AB'
            % params(1) = eta, params(2) = kappa
            % std_T1     = scale_factor(1) * sqrt(1 + params(1)^2 * params(2));
            % z_grid_std = linspace(-8 * std_T1, 8 * std_T1, N_grid)';
            [cdf_fT1, z_grid]    = lewis_FFT_digital(@cf_AB, M, dz, params, ...
                             scale_factor(1), 1, 'AB', 1);
            ft1        = simulate_from_cdf(cdf_fT1, z_grid, 1, N_sim);
    end

    % --- Forward price at t1 for each simulated path ---
    F_t1_T2 = forward + ft1;

    % --- Analytically price the inner call at t1 for each path ---
    % strikes are the shifted log-moneyness K2 - F_t1_T2
    strikes      = K2 - F_t1_T2;
    df_t1_t2     = discount_factors(2) / discount_factors(1);

    switch model
        case 'MA'
            call_price_t1 = df_t1_t2 * lewis_FFT_call(@cf_MA_FA, M, dz, params, ...
                                scale_factor, strikes, 1, 'MA');
        case 'GL'
            call_price_t1 = df_t1_t2 * lewis_FFT_call(@cf_increment_GL, M, dz, params, ...
                                scale_factor, strikes, 1, 'GL');
        case 'AB'
            call_price_t1 = df_t1_t2 * lewis_FFT_call(@cf_increment_AB, M, dz, params, ...
                                scale_factor, strikes, 1, 'AB');
    end
    
    put_price_t1   = call_price_t1 + df_t1_t2 * strikes;

    payoffs = max(call_price_t1, put_price_t1);
    discounted_payoffs = discount_factors(1)*payoffs;
    price   = mean(discounted_payoffs);
    std_err = std(discounted_payoffs) / sqrt(N_sim);
    CI      = [price - 1.96 * std_err, price + 1.96 * std_err];

end