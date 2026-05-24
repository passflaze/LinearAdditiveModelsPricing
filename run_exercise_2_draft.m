% data loading
clear,clc
addpath("Utilities/");
addpath("Functions/");
%%
callpath="Data/datacalls";
putpath="Data/dataputs";
expiryFile = "Data/Expiries_Futures.txt";
valueDate = datetime(2020,06,02);

[strikes, calls, puts, expiries] = readData(callpath, putpath, valueDate, expiryFile);

% bootstrap: synthetic discount factor and forward, one per maturity
nT = numel(expiries);
discount_factor = zeros(nT,1);
forward = zeros(nT,1);
R2 = zeros(nT,1);
for k = 1:nT
    [discount_factor(k), forward(k), R2(k)] = bootstrap(puts(k,:), calls(k,:), strikes);
end

% report
fprintf("\nBootstrap results - value date %s\n", string(valueDate, "yyyy-MM-dd"));
fprintf("%-12s  %10s  %10s  %8s\n", "Expiry", "D(t,T)", "F(t,T)", "R^2");
for k = 1:nT
    fprintf("%-12s  %10.6f  %10.4f  %8.4f\n", ...
        string(expiries(k), "yyyy-MM-dd"), ...
        discount_factor(k), forward(k), R2(k));
end

%% running ex 2

c_ATM = zeros(length(forward), 1); 

for i = 1:length(forward)
    % Extract data for current maturity
    current_calls   = calls(i, :);
    current_puts    = puts(i, :);
    
    c_ATM(i) = callATM(current_calls, current_puts, strikes, forward(i), discount_factor(i));
end

%% step 1 : calibrate sigma ATM
% since expiries are already bd adjusted we ignore the adjustment

yf = yearfrac(valueDate, expiries, 3);
sigma_ATM = sigmaATM(c_ATM, discount_factor, yf, expiries);
[moneyness_modified, c_mkt_calibration_normed, norm_factor] = moneyness_generator(forward,strikes,calls,puts,sigma_ATM,yf,discount_factor);
c_mkt_calibration = c_mkt_calibration_normed.*norm_factor;



%% MA calibration
% step 3 : calibrate market prices via fmincon (bounded, SSE objective)
% Symmetric starting point (gamma_MA = 1/alpha - 1/beta = 0) avoids initializing
% on a strongly skewed distribution. Bounds in [0.2, 10] keep 1/alpha and 1/beta
% in [0.1, 5], so gamma_MA stays within a sensible range and the closed-form
% bracket does not blow up numerically.
x0 = [1.0, 1.0];
lb = [0.2, 0.2];
ub = [10,  10];

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'interior-point', ...
    'FunctionTolerance',   1e-10, ...
    'StepTolerance',       1e-10, ...
    'OptimalityTolerance', 1e-8, ...
    'MaxFunctionEvaluations', 5000);

fprintf('\n--- Launching Volatility Surface Calibration via fmincon ---\n');
obj_fun_MA = @(x) objective_function_MA(x, discount_factor, yf, sigma_ATM, moneyness_modified, c_mkt_calibration);

[x_opt, fval, exitflag] = fmincon(obj_fun_MA, x0, [], [], [], [], lb, ub, [], options);

alpha_MA = x_opt(1);
beta_MA = x_opt(2);

fprintf('\nCalibration completed (exitflag = %d, SSE = %.6g).\n', exitflag, fval);
fprintf('Optimal Alpha: %.6f\n', alpha_MA);
fprintf('Optimal Beta : %.6f\n', beta_MA);
%% pricing trial
params = [alpha_MA,beta_MA];
price_MA_prova = price_MA(params,discount_factor,yf,sigma_ATM,moneyness_modified);

%% Calibration via GL
% step 3 : calibrate market prices via fmincon (bounded, SSE objective)
M = 15;
dz = 2.5e-3;
x0 = [0.9, 0.5];              
lb = [0.05, 0.05];             % positivity (replaces the discontinuous 1e10 penalty)
ub = [50,  50];                % loose upper cap

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'interior-point', ...
    'FunctionTolerance',   1e-10, ...
    'StepTolerance',       1e-10, ...
    'OptimalityTolerance', 1e-8, ...
    'MaxFunctionEvaluations', 5000);

fprintf('\n--- Launching Volatility Surface Calibration via fmincon ---\n');
obj_fun_GL = @(x) objective_function_GL(x, discount_factor, yf, sigma_ATM, moneyness_modified, c_mkt_calibration,M,dz);

[x_opt, fval, exitflag] = fmincon(obj_fun_GL, x0, [], [], [], [], lb, ub, [], options);

alpha_GL = x_opt(1);
beta_GL  = x_opt(2);

fprintf('\nCalibration completed (exitflag = %d, SSE = %.6g).\n', exitflag, fval);
fprintf('Optimal Alpha: %.6f\n', alpha_GL);
fprintf('Optimal Beta : %.6f\n', beta_GL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% diagnostica residui post-calibrazione
c_mod_MA = price_MA([alpha_MA, beta_MA], discount_factor, yf, sigma_ATM, moneyness_modified);
c_mod_GL = price_GL(alpha_GL, beta_GL, M, dz, discount_factor, sigma_ATM, yf, moneyness_modified);

res_MA = c_mkt_calibration - c_mod_MA;
res_GL = c_mkt_calibration - c_mod_GL;

% RMSE per maturity
fprintf('\n%-12s  %10s  %10s\n', 'Expiry', 'rMSE_MA', 'rMSE_GL');
for i = 1:nT
    r_MA_i = res_MA(i,:); r_MA_i = r_MA_i(~isnan(r_MA_i));
    r_GL_i = res_GL(i,:); r_GL_i = r_GL_i(~isnan(r_GL_i));
    fprintf('%-12s  %10.4f  %10.4f\n', string(expiries(i),'yyyy-MM-dd'), ...
        sqrt(mean(r_MA_i.^2)), sqrt(mean(r_GL_i.^2)));
end
