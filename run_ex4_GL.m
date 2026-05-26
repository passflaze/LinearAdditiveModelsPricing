% data loading
clear, clc
addpath("Utilities/");
addpath("Functions/");
addpath("ex3_GL/");
addpath("ex4_GL/"); 

callpath   = "Data/datacalls";
putpath    = "Data/dataputs";
expiryFile = "Data/Expiries_Futures.txt";
valueDate  = datetime(2020,06,02);

[strikes, calls, puts, expiries] = readData(callpath, putpath, valueDate, expiryFile);

% Bootstrap: discount factor e forward per ogni maturity
nT = numel(expiries);
discount_factor = zeros(nT,1);
forward         = zeros(nT,1);
R2              = zeros(nT,1);
for k = 1:nT
    [discount_factor(k), forward(k), R2(k)] = bootstrap(puts(k,:), calls(k,:), strikes);
end


%% PRICE CoC - GL model

alpha_GL = 0.40;
beta_GL  = 0.44;

% sigma ATM e normalizzazione GL
yf    = yearfrac(valueDate, expiries, 3);
c_ATM = zeros(length(forward), 1);
for i = 1:length(forward)
    c_ATM(i) = callATM(calls(i,:), puts(i,:), strikes, forward(i), discount_factor(i));
end
sigma_ATM      = sigmaATM(c_ATM, discount_factor, yf, expiries);
integrand_mean = @(x) pdf_GL(alpha_GL, beta_GL, x) .* x;
I0_opt         = sqrt(2*pi) * quadgk(integrand_mean, 0, inf);

% Tempi e sigma
t1 = yearfrac(valueDate, valueDate + calmonths(6), 3);
t2 = yearfrac(valueDate, valueDate + calyears(1),  3);
sigma_t1 = sigma_ATM(2) / I0_opt;
sigma_t2 = sigma_ATM(4) / I0_opt;

%DF 
B_t0_t1 = discount_factor(2);
B_t0_t2 = discount_factor(4);
B_t1_t2 = B_t0_t2/B_t0_t1; 

%prezzi call interna con FFT
[X_t1_grid, C_t1_grid] = call_GL_FFT(alpha_GL, beta_GL, t1, sigma_t1, t2, sigma_t2, B_t1_t2);

%simulo con MC tra t0 e t1 
[CDF_t1, zk_t1] = lewis_fft_cdf(@cf_increment_GL, alpha_GL, beta_GL, ...
                                0, 0, sigma_t1, t1);
[X1_sim, ~, ~]  = simulate_increments(zk_t1, CDF_t1, 1e5, 1);

%interpolare 
call_sim = interp1(X_t1_grid, C_t1_grid, X1_sim, 'pchip', 'extrap');

%payoff e price 
K1= B_t0_t1 * mean(call_sim); %CoC ATM rispetto call interna
payoff_CoC = max(call_sim - K1, 0);
price_CoC = B_t0_t1 * mean(payoff_CoC); 


%% estrapolazione smart 
p_minus_t1  = alpha_GL / (sigma_t1 * sqrt(t1));
p_plus_t1 = beta_GL  / (sigma_t1 * sqrt(t1));

[X1_sim_smart, ~, ~] = simulate_increments_smart(zk_t1, CDF_t1, 1e5, 1, p_minus_t1, p_plus_t1);
call_sim_smart = interp1(X_t1_grid, C_t1_grid, X1_sim_smart, 'pchip', 'extrap');
%payoff e price 
K1= B_t0_t1 * mean(call_sim_smart);  
payoff_CoC_smart = max(call_sim_smart - K1, 0);
price_CoC_smart = B_t0_t1 * mean(payoff_CoC_smart); 

%% PRICE PoP - GL model 
% put call parity su C_t1_grid 

P_t1_grid = C_t1_grid - B_t1_t2 * X_t1_grid; 
put_sim_smart = interp1(X_t1_grid, P_t1_grid, X1_sim_smart, 'pchip', 'extrap');

K1_put = B_t0_t1 * mean(put_sim_smart);
payoff_PoP_smart = max(K1_put - put_sim_smart, 0); 
price_PoP_smart = B_t0_t1 * mean(payoff_PoP_smart);


%% PRICE CHOOSER - GL model 

payoff_chooser_smart = max(call_sim_smart, put_sim_smart);
price_chooser_smart = B_t0_t1 * mean(payoff_chooser_smart);

%% print risultati
fprintf('\n--- RISULTATI FINALI ---\n');
fprintf('Prezzo CoC (ATM) : %.4f\n', price_CoC_smart);
fprintf('Prezzo PoP (ATM) : %.4f\n', price_PoP_smart);
fprintf('Prezzo Chooser   : %.4f\n', price_chooser_smart);


