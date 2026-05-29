function [price, IC] = price_chooser_GL_MC(T1, T2, alpha, beta, sigma_T1, sigma_T2, ...
                                           Nsim, forward, B_0_t2, N_grid)
% PRICE_CHOOSER_GL_MC  Monte Carlo price of a simple chooser under GL.
%
% At T1 the holder picks Call or Put on F with strike K = F(t0,T2) and
% maturity T2. Outer payoff at T1: max( C(T1; K, T2), P(T1; K, T2) ).
% Inner Call via Lewis-FFT (call_GL_FFT); inner Put via parity on the
% undiscounted quantities  Ptilde = Ctilde + x_money.
% Discounts collapse: B(0,T1)*B(T1,T2) = B(0,T2).


std_z     = sqrt(psi(1, alpha) + psi(1, beta));
std_T1    = sigma_T1 * sqrt(T1) * std_z;
x_grid_T1 = linspace(-10*std_T1, 10*std_T1, N_grid)';

cdf_fT1 = lewis_fft_cdf(alpha, beta, 0, T1, 0, sigma_T1, x_grid_T1, 1);
f_T1    = simulate_increments_smart(x_grid_T1, cdf_fT1, Nsim);
F_T1_T2 = forward + f_T1;

K       = forward;
x_money = K - F_T1_T2;
E_call  = call_GL_FFT(x_money, T1, T2, alpha, beta, sigma_T1, sigma_T2);
E_put   = E_call + x_money;

payoff  = max(E_call, E_put);            % B(T1,T2) factors out of max
price   = B_0_t2 * mean(payoff);         % B(0,T1)*B(T1,T2) = B(0,T2)
std_err = B_0_t2 * std(payoff) / sqrt(Nsim);
IC      = [price - 1.96*std_err, price + 1.96*std_err];

end