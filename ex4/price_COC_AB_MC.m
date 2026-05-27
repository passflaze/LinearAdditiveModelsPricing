function [price, IC] = price_COC_AB_MC(T1, T2, kAB, eta, sigma_T1, sigma_T2, ...
                                       Nsim, forward, B_0_t1, B_0_t2, N_grid, strike)
% PRICE_COC_AB_MC  Monte Carlo price of a Call-on-Call under Additive Bachelier.
%
% Outer payoff at T1:  max( C(T1; K2, T2) - K1, 0 )   with  K2 = F(t0, T2).
% Inner Call C(T1; K2, T2) is priced semi-analytically via Lewis-FFT on the
% increment CF phi_Delta = phi_{T2}/phi_{T1}, so the scheme is single-level
% MC (NOT nested): only f_{T1} is simulated; the inner expectation is done
% in closed-form by lewis_FFT_AB.
%
% Inputs
%   T1, T2   : outer / inner maturity, with T1 < T2
%   kAB, eta : AB parameters (constant in time, Prop.2.1)
%   sigma_T1 : AB scale at T1, = sigma_ATM(T1)/I_0   (Eq.15)
%   sigma_T2 : AB scale at T2, = sigma_ATM(T2)/I_0
%   Nsim     : number of MC paths
%   forward  : F(t0, T2)                              (also used as K2)
%   B_0_t1   : discount factor B(t0, T1)
%   B_0_t2   : discount factor B(t0, T2)
%   N_grid   : number of points for the f_{T1} marginal CDF grid
%   strike   : outer strike K1
%
% Outputs
%   price : MC estimate of the CoC price at t0
%   IC    : 95% confidence interval [lower, upper]

% 1) Simulate F(T1,T2) = F(t0,T2) + f_{T1,T2} via inverse-CDF on f_{T1,T1}
% Grid scaled to 8 stdev of f_{T1,T1}  (Lemma B.1: Var(z) = 1 + eta^2 * k).
% Lemma 2 of Forward.pdf: f_{T1,T2} = (B(0,T1)/B(0,T2)) * f_{T1,T1}.
std_T1    = sigma_T1 * sqrt(T1) * sqrt(1 + eta^2 * kAB);
x_grid_T1 = linspace(-8*std_T1, 8*std_T1, N_grid)';

cdf_fT1    = ccdf_AB_FFT(eta, kAB, 0, T1, 0, sigma_T1, x_grid_T1, 1);
f_T1       = sample_from_cdf(x_grid_T1, cdf_fT1, Nsim);
fwd_factor = B_0_t1 / B_0_t2;                          % = exp(int_{T1}^{T2} r_s ds)
F_T1_T2    = forward + fwd_factor * f_T1;              % Nsim x 1

%  2) Inner Call at T1 via Lewis-FFT on phi_{T2}/phi_{T1}(fwd_factor*.)
% lewis_FFT_AB returns the UNDISCOUNTED conditional expectation under Lemma 2;
% B(T1,T2) brings it back to a T1-price (deterministic rates).
K2      = forward;
x_money = K2 - F_T1_T2;                                % Nsim x 1
E_inner = lewis_FFT_AB(x_money, T1, T2, kAB, eta, sigma_T1, sigma_T2, fwd_factor);
B_T1_T2 = B_0_t2 / B_0_t1;
C_T1    = B_T1_T2 * E_inner;                           % Nsim x 1

% 3) Outer CoC payoff and MC estimator
payoff  = max(C_T1 - strike, 0);
price   = B_0_t1 * mean(payoff);
std_err = B_0_t1 * std(payoff) / sqrt(Nsim);
IC      = [price - 1.96*std_err, price + 1.96*std_err];

end
