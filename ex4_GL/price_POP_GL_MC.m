function [price, IC] = price_POP_GL_MC(T1, T2, alpha, beta, sigma_T1, sigma_T2, ...
                                       Nsim, forward, B_0_t1, B_0_t2, N_grid, strike)
% PRICE_POP_GL_MC  Monte Carlo price of a Put-on-Put under GL additive model.
%
% Outer payoff at T1:  max( K1 - P(T1; K2, T2), 0 )   with  K2 = F(t0, T2).
% Inner Put obtained from inner Call (call_GL_FFT) via put-call parity on
% undiscounted quantities:  Ptilde = Ctilde + x_money,  with
% x_money = K2 - F(T1,T2).

std_z     = sqrt(psi(1, alpha) + psi(1, beta));
std_T1    = sigma_T1 * sqrt(T1) * std_z;
x_grid_T1 = linspace(-10*std_T1, 10*std_T1, N_grid)';

cdf_fT1 = lewis_fft_cdf(alpha, beta, 0, T1, 0, sigma_T1, x_grid_T1, 1);
f_T1    = simulate_increments_smart(x_grid_T1, cdf_fT1, Nsim);
F_T1_T2 = forward + f_T1;

K2      = forward;
x_money = K2 - F_T1_T2;
E_call  = call_GL_FFT(x_money, T1, T2, alpha, beta, sigma_T1, sigma_T2);
E_put   = E_call + x_money;                              % parity (undiscounted)

B_T1_T2 = B_0_t2 / B_0_t1;
P_T1    = B_T1_T2 * E_put;                               % T1-price of inner Put

payoff  = max(strike - P_T1, 0);
price   = B_0_t1 * mean(payoff);
std_err = B_0_t1 * std(payoff) / sqrt(Nsim);
IC      = [price - 1.96*std_err, price + 1.96*std_err];

end