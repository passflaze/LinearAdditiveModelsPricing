% run_ex4_GL.m
% Esercizio 4: CoC, PoP, Chooser sotto modello Generalized Logistic (GL).

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

%% Bootstrap: discount factor, forward, sigma_ATM per ogni maturity
nT = numel(expiries);
discount_factor = zeros(nT,1);
forward         = zeros(nT,1);
R2              = zeros(nT,1);
call_atm        = zeros(nT,1);
for k = 1:nT
    [discount_factor(k), forward(k), R2(k)] = bootstrap(puts(k,:), calls(k,:), strikes);
    call_atm(k) = callATM(calls(k,:), puts(k,:), strikes, forward(k), discount_factor(k));
end

act_365   = 3;
yf        = yearfrac(valueDate, expiries, act_365);
sigma_atm = sigmaATM(call_atm, discount_factor, yf, expiries);

%% Parametri GL (da calibrazione esercizio 2) e scala sigma_t = sigma_ATM/I_0
alpha_GL = 0.40;
beta_GL  = 0.44;

integrand_mean = @(x) pdf_GL(alpha_GL, beta_GL, x) .* x;
I0_opt         = sqrt(2*pi) * quadgk(integrand_mean, 0, inf);

sigma_t = sigma_atm / I0_opt;          % sigma_t(iT) = sigma_ATM(iT)/I_0

%% Selezione maturity T1 = 6m, T2 = 1y 
T1       = yf(2);       T2       = yf(4);
sigma_T1 = sigma_t(2);  sigma_T2 = sigma_t(4);

seed   = 1234;
Nsim   = 1e6;
N_grid = 300;

%% --- 4a Call-on-Call (CoC) ---------------------------------------------
K1 = 1;                                  

rng(seed);
[CoC_price, CoC_IC] = price_COC_GL_MC(T1, T2, alpha_GL, beta_GL, ...
                                      sigma_T1, sigma_T2, Nsim, ...
                                      forward(4), discount_factor(2), ...
                                      discount_factor(2), N_grid, K1);

% Sanity: K1=0 collassa l'outer max -> vanilla ATM-forward Call at T2.
rng(seed);
[CoC_K0, ~] = price_COC_GL_MC(T1, T2, alpha_GL, beta_GL, ...
                              sigma_T1, sigma_T2, Nsim, ...
                              forward(4), discount_factor(2), ...
                              discount_factor(2), N_grid, 0);
G0       = G0_GL(alpha_GL, beta_GL);                    % E[z^+] normalizzato
C_T2_ATM = discount_factor(4) * sigma_atm(4) * sqrt(T2) * G0;

fprintf('\n--- 4a CoC GL (K1=%.2f, K2=F(t0,T2)=%.4f, T1=%.2fy, T2=%.2fy) ---\n', ...
        K1, forward(iT2), T1, T2);
fprintf('MC price       : %.4f   95%% CI [%.4f, %.4f]\n', CoC_price, CoC_IC(1), CoC_IC(2));
fprintf('Sanity (K1=0)  : MC %.4f  vs  FFT Call_T2_ATM %.4f   (rel err %.2e)\n', ...
        CoC_K0, C_T2_ATM, abs(CoC_K0 - C_T2_ATM)/C_T2_ATM);

%% --- 4a Put-on-Put (PoP) -----------------------------------------------
rng(seed);
[PoP_price, PoP_IC] = price_POP_GL_MC(T1, T2, alpha_GL, beta_GL, ...
                                      sigma_T1, sigma_T2, Nsim, ...
                                      forward(4), discount_factor(2), ...
                                      discount_factor(4), N_grid, K1);

fprintf('\n--- 4a PoP GL (K1=%.2f, K2=F(t0,T2)=%.4f) ---\n', K1, forward(iT2));
fprintf('MC price       : %.4f   95%% CI [%.4f, %.4f]\n', PoP_price, PoP_IC(1), PoP_IC(2));

%% --- 4b Chooser --------------------------------------------------------
rng(seed);
[ch_price, ch_IC] = price_chooser_GL_MC(T1, T2, alpha_GL, beta_GL, ...
                                        sigma_T1, sigma_T2, Nsim, ...
                                        forward(iT2), discount_factor(iT2), N_grid);
 
% Sanity (Stochastic-Reset trick, come collega):
% Chooser = B(0,T2) * G(0) * [sigma_ATM(T2)*sqrt(T2) + sigma_ATM(T1)*sqrt(T1)]
ch_analytic = discount_factor(iT2) * G0 * ...
              ( sigma_atm(iT2)*sqrt(T2) + sigma_atm(iT1)*sqrt(T1) );
 
fprintf('\n--- 4b Chooser GL (K=F(t0,T2)=%.4f, T1=%.2fy, T2=%.2fy) ---\n', ...
        forward(iT2), T1, T2);
fprintf('MC price       : %.4f   95%% CI [%.4f, %.4f]\n', ch_price, ch_IC(1), ch_IC(2));
fprintf('Analytic (SR)  : %.4f                            (rel err %.2e)\n', ...
        ch_analytic, abs(ch_price - ch_analytic)/ch_analytic);


%% ----------------------------------------------------------------------
function g0 = G0_GL(alpha, beta)
% G(0) = E[ z^+ ] for z ~ GL(alpha,beta) (standardized via cf_GL).
% Used in the analytical sanity checks (Stochastic-Reset trick).
    integrand = @(x) pdf_GL(alpha, beta, x) .* x;
    g0 = quadgk(integrand, 0, inf);
end