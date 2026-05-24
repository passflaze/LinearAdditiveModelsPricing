function [price, IC]= price_AB_MC(T1, T2, kAB, eta, sigma_T1, sigma_T2, Nsim, forward, discount_factor,x_grid, strike)

cdf_1= ccdf_AB_FFT(eta, kAB, 0, T1, 0, sigma_T1, x_grid);
Z_1 = sample_from_cdf(x_grid, cdf_1, Nsim);

% F(t1,t2)= forward+xt1
f1_2= forward * ones(Nsim, 1) + Z_1;

ccdf= ccdf_AB_FFT(eta, kAB, T1, T2, sigma_T1, sigma_T2, x_grid);
Z_2_1 = sample_from_cdf(x_grid, ccdf, Nsim);

% St2=F(t2,t2)= F(t0,t2)+xt1+deltax

f2_2= forward + Z_1 + Z_2_1;

payoff = max(f2_2 - strike*f1_2, 0);
price = discount_factor * mean(payoff);
std_err = discount_factor * std(payoff) / sqrt(Nsim);
IC = [price - 1.96*std_err, price + 1.96*std_err];


end