function [price, IC, diag] = price_fwd_start_MC(model, params, T1, T2, ...
        sigma_T1, sigma_T2, Nsim, forward, B_0_t1, B_0_t2, x_grid, strike, seed)
% PRICE_FWD_START_MC  Monte Carlo price of the forward-start option
%   payoff = [ S(T2) - strike * F(T1,T2) ]_+
% under a Linear Additive model ('AB' or 'GL'). Model-agnostic generalisation
% of ex3_utilities/price_AB_MC.
%
% Strategy (exploits additivity x_T2 = x_T1 + (x_T2 - x_T1), independent
% increments) and the Lemma 2 forward rescaling (Forward.pdf):
%   1. sample Z1 = f_{T1,T1}            from the marginal CDF   (0  -> T1)
%   2. sample W  = f_{T2,T2}-f_{T1,T2}  from the conditional CDF (T1 -> T2)
%   3. f_{T1,T2} = fwd_factor * f_{T1,T1},  fwd_factor = B(0,T1)/B(0,T2):
%        F(T1,T2) = forward + fwd_factor * Z1
%        S(T2)    = forward + fwd_factor * Z1 + W
% Both CDFs go through the SAME pipeline: ccdf_increment_FFT (Lewis two-shift
% + Lemma 2) -> tail_adjustment (clean/refine/exp-tails) -> simulate_from_cdf
% (monotone PCHIP inversion).
%
% INPUTS
%   model            : 'AB' or 'GL'
%   params           : model parameter struct (see model_marginal_cf)
%   T1, T2           : reset and maturity (year fractions)
%   sigma_T1,sigma_T2: calibrated scales at T1, T2 (= sigma_ATM/I_0)
%   Nsim             : number of MC paths
%   forward          : F(0, T2)
%   B_0_t1           : B(0, T1)   (Lemma 2 forward rescaling)
%   B_0_t2           : B(0, T2)   (payoff received at T2)
%   x_grid           : evaluation grid for the CDFs ($-increments, column)
%   strike           : K (the K2 of the forward-start option)
%   seed             : (optional) RNG seed for reproducibility
%
% OUTPUTS
%   price : B(0,T2) * E[(S(T2) - K*F(T1,T2))_+]
%   IC    : 95% confidence interval [lo, hi]
%   diag  : struct with the cleaned conditional CDF (x_cond, cdf_cond) for
%           plotting / FFT-reference cross-checks.

    if nargin >= 13 && ~isempty(seed)
        rng(seed);
    end

    x_grid     = x_grid(:);
    fwd_factor = B_0_t1 / B_0_t2;          % Lemma 2: = exp(int_{T1}^{T2} r_s ds)

    % --- marginal CDF (0 -> T1) -------------------------------------------
    cdf_1          = ccdf_increment_FFT(model, params, 0, T1, 0, sigma_T1, x_grid, 1);
    [cdf_1c, x_1c] = tail_adjustment(x_grid, cdf_1, 10);
    Z1             = simulate_from_cdf(cdf_1c, x_1c, true, Nsim);

    % --- conditional CDF (T1 -> T2) under Lemma 2 -------------------------
    cdf_2          = ccdf_increment_FFT(model, params, T1, T2, sigma_T1, sigma_T2, x_grid, fwd_factor);
    [cdf_2c, x_2c] = tail_adjustment(x_grid, cdf_2, 10);
    W              = simulate_from_cdf(cdf_2c, x_2c, true, Nsim);

    % --- Lemma 2 reconstruction of forward / spot -------------------------
    F_T1_T2 = forward + fwd_factor * Z1;          % F(T1, T2)
    S_T2    = forward + fwd_factor * Z1 + W;      % S(T2) = F(T2, T2)

    % --- payoff and MC estimator ------------------------------------------
    payoff  = max(S_T2 - strike * F_T1_T2, 0);
    price   = B_0_t2 * mean(payoff);
    std_err = B_0_t2 * std(payoff) / sqrt(Nsim);
    IC      = [price - 1.96*std_err, price + 1.96*std_err];

    diag = struct('x_cond', x_2c, 'cdf_cond', cdf_2c, 'W', W);
end
