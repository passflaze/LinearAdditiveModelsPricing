function [price, IC] = price_AB_MC(T1, T2, kAB, eta, sigma_T1, sigma_T2, ...
                                    Nsim, forward, discount_factor, x_grid, strike)
% PRICE_AB_MC  MC price of the forward start option [S(T2) - K*F(T1,T2)]_+
% under the Additive Bachelier model.
%
% Strategy (exploits additivity of {x_t}: x_T2 = x_T1 + (x_T2-x_T1) with
% independent increments):
%   1. sample Z_1   ~ x_{T1}              from the marginal CDF (0 -> T1)
%   2. sample Z_2_1 ~ x_{T2} - x_{T1}     from the conditional CDF (T1 -> T2)
%   3. linear-forward relation (Forward.pdf, Lemma 2): F(t, T2) = F(0,T2) + x_t
%        F(T1, T2) = forward + Z_1
%        S(T2)     = F(T2,T2) = forward + Z_1 + Z_2_1
%
% Inputs
%   T1, T2          : reset and maturity (year fractions)
%   kAB, eta        : AB model parameters
%   sigma_T1,sigma_T2: calibrated sigma_t at T1 and T2
%   Nsim            : number of MC paths
%   forward         : F(0, T2)   - initial forward at T2 
%   discount_factor : B(0, T2)   - payoff is received at T2
%   x_grid          : grid for the conditional CDF T1 -> T2
%   strike          : K (the K_2 of the forward start option)
%
% Outputs
%   price : MC estimator of B(0,T2) * E_0[(S(T2) - K*F(T1,T2))_+]
%   IC    : 95% confidence interval [lo, hi]

% --- marginal CDF (0 -> T1) for the T2-forward ----------------------------
% Lemma 2 (Forward.pdf, Sec. 1.2): x_{T1}^{T2} = D(T1,T2) * x_{T1}^{T1}
% where D = B(0,T1)/B(0,T2). Equivalently, x_{T1}^{T2} ~ AB(eta,k,sigma_T2,T1)
% so we use sigma_T2 here to get Z_1 ~ x_{T1}^{T2} directly.
% AB variance: Var(x_T) = sigma_T^2 * T * (1 + eta^2 * k).
std_T1    = sigma_T2 * sqrt(T1) * sqrt(1 + eta^2 * kAB);
x_grid_T1 = linspace(-8*std_T1, 8*std_T1, numel(x_grid))';

cdf_1 = ccdf_AB_FFT(eta, kAB, 0, T1, 0, sigma_T2, x_grid_T1);
Z_1   = sample_from_cdf(x_grid_T1, cdf_1, Nsim);

% --- conditional CDF (T1 -> T2) -------------------------------------------
ccdf  = ccdf_AB_FFT(eta, kAB, T1, T2, sigma_T1, sigma_T2, x_grid);
Z_2_1 = sample_from_cdf(x_grid, ccdf, Nsim);

% --- reconstruct forward / spot via linear-forward relation ---------------
f1_2 = forward + Z_1;            % F(T1, T2)
f2_2 = forward + Z_1 + Z_2_1;    % S(T2) = F(T2, T2)

% --- payoff and MC estimator ----------------------------------------------
payoff  = max(f2_2 - strike * f1_2, 0);
price   = discount_factor * mean(payoff);
std_err = discount_factor * std(payoff) / sqrt(Nsim);
IC      = [price - 1.96*std_err, price + 1.96*std_err];

end
