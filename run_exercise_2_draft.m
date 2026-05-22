% data loading
clear,clc
addpath("Utilities/");
addpath("Functions/");

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
    [discount_factor(k), forward(k), R2(k)] = bootstrap_forward(puts(k,:), calls(k,:), strikes);
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

%% step 2 : defining pricefunctions 
 
M = length(forward);
N = length(strikes);

moneyness_modified = NaN(M, N);
c_mkt_calibration  = NaN(M, N);

for i = 1:M
    mask = strikes > forward(i);

    moneyness_modified(i, mask) = (strikes(mask) - forward(i)) / (sigma_ATM(i) * sqrt(yf(i)));
   
    c_mkt_calibration(i, mask) = calls(i, mask);
end

%%
% step 3 : calibrate market prices 
x0 = [0.5,0.5];
% 4. Set optimization options for fminsearch
options = optimset('Display', 'iter', 'TolX', 1e-6, 'TolFun', 1e-6);

% 5. Run fminsearch (Notice: no bounds lb, ub are passed here!)
fprintf('\n--- Launching Volatility Surface Calibration via fminsearch ---\n');
obj_fun = @(x) objective_function_MA(x, discount_factor, yf, sigma_ATM, moneyness_modified, c_mkt_calibration);

% Run the simplex optimizer
x_opt = fminsearch(obj_fun, x0, options);

% 6. Map optimal parameters back into the structured model format
alpha_MA = x_opt(1);
beta_MA  = x_opt(2);

fprintf('\nCalibration completed.\n');
fprintf('Optimal Alpha: %.6f\n', alpha_MA);
fprintf('Optimal Beta : %.6f\n', beta_MA);


%% GL model : calibration 

% step 3 : calibrate market prices 
x0_GL = [0.5,0.5];
M = 15;
dz = 2.5e-3;
% 4. Set optimization options for fminsearch
options = optimset('Display', 'iter', 'TolX', 1e-6, 'TolFun', 1e-6);

% 5. Run fminsearch (Notice: no bounds lb, ub are passed here!)
fprintf('\n--- Launching Volatility Surface Calibration via fminsearch ---\n');
obj_fun_GL = @(x) objective_function_GL(x, discount_factor, yf, sigma_ATM, moneyness_modified, c_mkt_calibration,M,dz);

% Run the simplex optimizer
x_opt_GL = fminsearch(obj_fun_GL, x0_GL, options);

% 6. Map optimal parameters back into the structured model format
alpha_GL = x_opt_GL(1);
beta_GL  = x_opt_GL(2);

fprintf('\nCalibration completed.\n');
fprintf('Optimal Alpha: %.6f\n', alpha_GL);
fprintf('Optimal Beta : %.6f\n', beta_GL);

%% prova fourier pricing (funziona)
alpha_GL = 0.9;
beta_GL = 1.5;
M = 15;
dz = 2.5e-3;

prova = price_GL(alpha_GL, beta_GL, M, dz, discount_factor, sigma_ATM, yf, moneyness_modified);

%% dà problemi !!
alpha_GL = 0.5;
beta_GL = 0.5;
shift = beta_GL / 2;
integrand = @(t) (1) ./ ((t - 1i*shift).^2) .* cf_GL(t - 1i*shift, alpha_GL, beta_GL);
I0 = (-1 / sqrt(2*pi)) * real(quadgk(integrand, -inf, inf));
