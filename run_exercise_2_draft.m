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
x0 = [0.5, 0.5];               
lb = [1e-4, 1e-4];             % positivity (replaces the discontinuous 1e10 penalty)
ub = [50,  50];                % loose upper cap

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'interior-point', ...
    'FunctionTolerance',   1e-10, ...
    'StepTolerance',       1e-10, ...
    'OptimalityTolerance', 1e-8, ...
    'MaxFunctionEvaluations', 5000);

fprintf('\n--- Launching Volatility Surface Calibration via fmincon ---\n');
obj_fun_MA = @(x) objective_function_MA(x, discount_factor, yf, sigma_ATM, moneyness_modified, c_mkt_calibration_normed);

[x_opt, fval, exitflag] = fmincon(obj_fun_MA, x0, [], [], [], [], lb, ub, [], options);

alpha_MA = x_opt(1);
beta_MA = x_opt(2);

fprintf('\nCalibration completed (exitflag = %d, SSE = %.6g).\n', exitflag, fval);
fprintf('Optimal Alpha: %.6f\n', alpha_MA);
fprintf('Optimal Beta : %.6f\n', beta_MA);
%% pricing trial
params = [0.001295,32.012298];
price_MA_prova = price_MA(params,discount_factor,yf,sigma_ATM,moneyness_modified);

%% Calibration via GL
% step 3 : calibrate market prices via fmincon (bounded, SSE objective)
M = 15;
dz = 2.5e-3;
x0 = [0.5, 0.5];               % asymmetric guess: WTI Covid period -> left tail heavier
lb = [1e-4, 1e-4];             % positivity (replaces the discontinuous 1e10 penalty)
ub = [50,  50];                % loose upper cap

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'interior-point', ...
    'FunctionTolerance',   1e-10, ...
    'StepTolerance',       1e-10, ...
    'OptimalityTolerance', 1e-8, ...
    'MaxFunctionEvaluations', 5000);

fprintf('\n--- Launching Volatility Surface Calibration via fmincon ---\n');
obj_fun_GL = @(x) objective_function_GL(x, discount_factor, yf, sigma_ATM, moneyness_modified, c_mkt_calibration_normed,M,dz);

[x_opt, fval, exitflag] = fmincon(obj_fun_GL, x0, [], [], [], [], lb, ub, [], options);

alpha_GL = x_opt(1);
beta_GL  = x_opt(2);

fprintf('\nCalibration completed (exitflag = %d, SSE = %.6g).\n', exitflag, fval);
fprintf('Optimal Alpha: %.6f\n', alpha_GL);
fprintf('Optimal Beta : %.6f\n', beta_GL);

%% prova output

provaGLprezzo = price_GL(alpha_GL, beta_GL, M, dz, discount_factor, sigma_ATM, yf, moneyness_modified);
%%
prova = cf_GL(xgrid-1i*beta_GL,alpha_GL, beta_GL);