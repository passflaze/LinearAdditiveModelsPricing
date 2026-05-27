function [price, IC] = price_AB_MC(T1, T2, kAB, eta, sigma_T1, sigma_T2, ...
                                    Nsim, forward, B_0_t1, B_0_t2, x_grid, strike)
% PRICE_AB_MC  MC price of the forward start option [S(T2) - K*F(T1,T2)]_+
% under the Additive Bachelier model.
%
% Strategy (exploits additivity of {x_t}: x_T2 = x_T1 + (x_T2-x_T1) with
% independent increments):
%   1. sample Z_1   ~ f_{T1,T1}           from the marginal CDF (0 -> T1)
%   2. sample Z_2_1 ~ f_{T2,T2}-f_{T1,T1} from the conditional CDF (T1 -> T2)
%   3. linear-forward relation (Forward.pdf, Lemma 2):
%        F(T1, T2) = F(0,T2) + f_{T1,T2}, with f_{T1,T2} = (B(0,T1)/B(0,T2))*f_{T1,T1}
%                  = forward + (B_0_t1/B_0_t2) * Z_1
%        S(T2)     = F(T2,T2) = F(0,T2) + f_{T2,T2}
%                  = forward + Z_1 + Z_2_1
%      The B_0_t1/B_0_t2 factor on Z_1 corrects for the fact that the linear
%      instantaneous forward depends on the expiry (Lemma 2): the AB process
%      Z_1 is calibrated to f_{T1,T1}, so f_{T1,T2} differs by exp(int_{T1}^{T2} r_s ds).
%
% Inputs
%   T1, T2          : reset and maturity (year fractions)
%   kAB, eta        : AB model parameters
%   sigma_T1,sigma_T2: calibrated sigma_t at T1 and T2
%   Nsim            : number of MC paths
%   forward         : F(0, T2)   - initial forward at T2
%   B_0_t1          : B(0, T1)   - discount factor to T1 (Lemma 2 forward rescaling)
%   B_0_t2          : B(0, T2)   - payoff is received at T2
%   x_grid          : grid for the conditional CDF T1 -> T2
%   strike          : K (the K_2 of the forward start option)
%
% Outputs
%   price : MC estimator of B(0,T2) * E_0[(S(T2) - K*F(T1,T2))_+]
%   IC    : 95% confidence interval [lo, hi]

% --- marginal CDF (0 -> T1) -----------------------------------------------
% The AB model defines a single process {f_t} with time-varying vol sigma_s
% (Baviera-Massaria 2026, Eq. 4). F(t,T) = F(0,T) + f_t for any T, so f_{T1}
% has the calibrated marginal with sigma_T1 — no cross-maturity rescaling needed.
% AB variance: Var(f_T) = sigma_T^2 * T * (1 + eta^2 * k).
std_T1    = sigma_T1 * sqrt(T1) * sqrt(1 + eta^2 * kAB);
x_grid_T1 = linspace(-8*std_T1, 8*std_T1, numel(x_grid))';

cdf_1 = ccdf_AB_FFT(eta, kAB, 0, T1, 0, sigma_T1, x_grid_T1,0);
Z_1   = sample_from_cdf(x_grid_T1, cdf_1, Nsim);

% --- conditional CDF (T1 -> T2) -------------------------------------------
ccdf  = ccdf_AB_FFT(eta, kAB, T1, T2, sigma_T1, sigma_T2, x_grid, 0);
Z_2_1 = sample_from_cdf(x_grid, ccdf, Nsim);


% --- reconstruct forward / spot via linear-forward relation ---------------
% Lemma 2 of Forward.pdf: f_{T1,T2} = (B(0,T1)/B(0,T2)) * f_{T1,T1}.
% Z_1 is calibrated to f_{T1,T1}, so the T2-forward at T1 uses the rate factor.
fwd_factor = B_0_t1 / B_0_t2;     % = exp(int_{T1}^{T2} r_s ds)
f1_2 = forward + fwd_factor * Z_1;   % F(T1, T2)
f2_2 = forward + Z_1 + Z_2_1;        % S(T2) = F(T2, T2)

% --- payoff and MC estimator ----------------------------------------------
payoff  = max(f2_2 - strike * f1_2, 0);
price   = B_0_t2 * mean(payoff);
std_err = B_0_t2 * std(payoff) / sqrt(Nsim);
IC      = [price - 1.96*std_err, price + 1.96*std_err];

end
