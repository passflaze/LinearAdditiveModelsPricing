function [price, IC] = price_chooser_AB_MC(T1, T2, kAB, eta, sigma_T1, sigma_T2, ...
                                           Nsim, forward, B_0_t2, x_grid)
% PRICE_CHOOSER_AB_MC  Monte Carlo price of a simple chooser under Additive Bachelier.
%
% At T1 the holder picks Call or Put on F with strike K = F(t0,T2) and maturity T2.
% Outer payoff at T1: max( C(T1; K, T2), P(T1; K, T2) ).
% Lewis-FFT gives the inner Call from phi_{T2}/phi_{T1}; the inner Put follows
% from put-call parity on undiscounted quantities (Ptilde = Ctilde + x).
% Discounts collapse: B(0,T1)*B(T1,T2) = B(0,T2), so only B(0,T2) appears.
%
% Inputs
%   T1, T2     : choice date / final maturity (T1 < T2)
%   kAB, eta   : AB parameters
%   sigma_T1   : AB scale at T1, = sigma_ATM(T1)/I_0
%   sigma_T2   : AB scale at T2, = sigma_ATM(T2)/I_0
%   Nsim       : number of MC paths
%   forward    : F(t0, T2)                       (also used as strike K)
%   B_0_t2     : discount factor B(t0, T2)
%   x_grid     : moneyness grid for the inverse-CDF sampler ($ column)
%
% Outputs
%   price : MC estimate of the chooser price at t0
%   IC    : 95% confidence interval [lower, upper]

% Simulate F(T1,T2) via marginal CDF of f_{T1} (T1==0 branch of ccdf_AB_FFT)
% Adjust the grid to the marginal variance of T1 to avoid tail mass warnings
std_T1    = sigma_T1 * sqrt(T1) * sqrt(1 + eta^2 * kAB);
x_grid_T1 = linspace(-8*std_T1, 8*std_T1, numel(x_grid))';

cdf_fT1 = ccdf_AB_FFT(eta, kAB, 0, T1, 0, sigma_T1, x_grid_T1, 1);
f_T1    = sample_from_cdf(x_grid_T1, cdf_fT1, Nsim);
F_T1_T2 = forward + f_T1;                              % Nsim x 1

% Inner Call (undiscounted) via Lewis on phi_{T2}/phi_{T1}
K       = forward;                                     % ATM-forward strike
x_money = K - F_T1_T2;                                 % Nsim x 1
E_call  = lewis_FFT_AB(x_money, T1, T2, kAB, eta, sigma_T1, sigma_T2);

% Inner Put (undiscounted) via put-call parity
E_put   = E_call + x_money;

% B(T1,T2) factors out of max; combined with B(0,T1) gives B(0,T2)
payoff  = max(E_call, E_put);
price   = B_0_t2 * mean(payoff);
std_err = B_0_t2 * std(payoff) / sqrt(Nsim);
IC      = [price - 1.96*std_err, price + 1.96*std_err];

end
