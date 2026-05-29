function [price, IC] = price_COC_GL_MC(T1, T2, alpha, beta, sigma_T1, sigma_T2, ...
                                       Nsim, forward, B_0_t1, B_0_t2, N_grid, strike)
% PRICE_COC_GL_MC  Monte Carlo price of a Call-on-Call under GL additive model.
%
% Outer payoff at T1:  max( C(T1; K2, T2) - K1, 0 )   with  K2 = F(t0, T2).
% Inner Call C(T1; K2, T2) is priced semi-analytically via Lewis-FFT
% (call_GL_FFT) on the increment CF phi_Delta = phi_{T2}/phi_{T1}, so the
% scheme is single-level MC (NOT nested): only f_{T1} is simulated; the
% inner expectation is closed-form.
%
% Inputs
%   T1, T2         : outer / inner maturity, with T1 < T2
%   alpha, beta    : GL shape parameters
%   sigma_T1       : GL scale at T1, = sigma_ATM(T1)/I_0
%   sigma_T2       : GL scale at T2, = sigma_ATM(T2)/I_0
%   Nsim           : number of MC paths
%   forward        : F(t0, T2)                          (also used as K2)
%   B_0_t1, B_0_t2 : discount factors
%   N_grid         : points for the f_{T1} marginal CDF grid
%   strike         : outer strike K1
%
% Outputs
%   price : MC estimate of the CoC price at t0
%   IC    : 95% confidence interval

% 1) Simulate F(T1,T2) = F(t0,T2) + f_{T1} via inverse-CDF
% GL is heavy-tailed: widen the grid. Variance of normalized z is
% psi'(alpha) + psi'(beta) (trigamma).
std_z     = sqrt(psi(1, alpha) + psi(1, beta));
std_T1    = sigma_T1 * sqrt(T1) * std_z;
x_grid_T1 = linspace(-10*std_T1, 10*std_T1, N_grid)';

cdf_fT1 = lewis_fft_cdf(alpha, beta, 0, T1, 0, sigma_T1, x_grid_T1, 1);
f_T1    = simulate_increments_smart(x_grid_T1, cdf_fT1, Nsim);
F_T1_T2 = forward + f_T1;                                % Nsim x 1

% 2) Inner Call at T1 via Lewis-FFT on phi_{T2}/phi_{T1}
K2      = forward;
x_money = K2 - F_T1_T2;                                  % Nsim x 1
E_inner = call_GL_FFT(x_money, T1, T2, alpha, beta, sigma_T1, sigma_T2);
B_T1_T2 = B_0_t2 / B_0_t1;
C_T1    = B_T1_T2 * E_inner;                             % Nsim x 1 (T1-price)

% 3) Outer CoC payoff and MC estimator
payoff  = max(C_T1 - strike, 0);
price   = B_0_t1 * mean(payoff);
std_err = B_0_t1 * std(payoff) / sqrt(Nsim);
IC      = [price - 1.96*std_err, price + 1.96*std_err];

end