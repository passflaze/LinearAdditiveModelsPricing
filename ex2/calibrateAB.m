function [k,eta,sigma_t,rMSE] = calibrateAB(marketNormMoneyness, marketNormCalls, sigmaATM)

% CALIBRATEAB  Calibrate the Additive Bachelier model parameters (k, eta) to market data.

c_mod= @(k,eta) call_AB_FFT(marketNormMoneyness, k, eta);

% Objective function: mean squared error between model and market normalized call prices
objective = @(params) sum(abs((c_mod(params(1), params(2)) - marketNormCalls).^2));     

% Initial guess for [k, eta]
initial_guess = [1.0, 0.2];

% Optimization: k > 0 (AB model requires k in R+), eta free in R
lb = [1e-6, -Inf];
ub = [ Inf,  Inf];
options = optimoptions('fmincon', 'Display', 'iter', ...
                       'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-10);
% Run optimization
[param_opt, MSE] = fmincon(objective, initial_guess, [], [], [], [], ...
                           lb, ub, [], options);

k = param_opt(1);
eta = param_opt(2); 
I_0=I0(0,k,eta);
sigma_t = sigmaATM /I_0;

rMSE = sqrt(MSE/numel(marketNormCalls)); % root mean squared error

end






















