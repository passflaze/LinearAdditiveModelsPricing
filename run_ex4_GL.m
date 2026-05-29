% run_ex4_GL.m
% Esercizio 4: CoC, PoP, Chooser sotto modello Generalized Logistic (GL).
% Schema single-level: si simula solo f_{T1}; la call/put interna a T1 e'
% calcolata in forma semi-chiusa via Lewis-FFT sulla CF dell'incremento
% phi_Delta = phi_{T2}/phi_{T1}.

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

%% Parametri GL (da calibrazione esercizio 2)
alpha_GL = 0.40;
beta_GL  = 0.44;

% Normalizzazione I_0 e scala normalizzata sigma_t = sigma_ATM / I_0.
% I_0 = sqrt(2*pi) * E[z^+] = sqrt(2*pi) * G0_GL.
integrand_mean = @(x) pdf_GL(alpha_GL, beta_GL, x) .* x;
I0_opt         = sqrt(2*pi) * quadgk(integrand_mean, 0, inf);
sigma_t        = sigma_atm / I0_opt;          % vettore: sigma_t(iT) = sigma_ATM(iT)/I_0

% G0 = E[z^+] = prezzo normalizzato call ATM-forward (numero puro)
G0 = G0_GL(alpha_GL, beta_GL);

%% Selezione maturity: T1 = 6m, T2 = 1y -> indici 2 e 4 nel TUO dataset
iT1 = 2;  iT2 = 4;
T1       = yf(iT1);       T2       = yf(iT2);
sigma_T1 = sigma_t(iT1);  sigma_T2 = sigma_t(iT2);

seed   = 1234;
Nsim   = 1e6;
N_grid = 300;

%% ====================== 4a Call-on-Call (CoC) =========================
% Prezzo VERO del contratto: strike esterno K1 fissato (placeholder = 1).
K1 = 1;

rng(seed);
[CoC_price, CoC_IC] = price_COC_GL_MC(T1, T2, alpha_GL, beta_GL, ...
                                      sigma_T1, sigma_T2, Nsim, ...
                                      forward(iT2), discount_factor(iT1), ...
                                      discount_factor(iT2), N_grid, K1);

% --- Sanity check (calcolo SEPARATO, K1 = 0) --------------------------
% Con K1=0 il max esterno collassa: CoC -> vanilla call ATM-forward a T2.
% Riferimento analitico (USA sigma_t, NON sigma_atm: la scala dei pricer
% e' sigma_t = sigma_ATM/I0; mischiarle introduce un fattore I0 ~ 3.38).
rng(seed);
[CoC_K0, ~] = price_COC_GL_MC(T1, T2, alpha_GL, beta_GL, ...
                              sigma_T1, sigma_T2, Nsim, ...
                              forward(iT2), discount_factor(iT1), ...
                              discount_factor(iT2), N_grid, 0);
C_T2_ATM = discount_factor(iT2) * sigma_t(iT2) * sqrt(T2) * G0;

fprintf('\n--- 4a CoC GL (K1=%.2f, K2=F(t0,T2)=%.4f, T1=%.2fy, T2=%.2fy) ---\n', ...
        K1, forward(iT2), T1, T2);
fprintf('MC price       : %.4f   95%% CI [%.4f, %.4f]\n', CoC_price, CoC_IC(1), CoC_IC(2));
fprintf('Sanity (K1=0)  : MC %.4f  vs  Call_T2_ATM %.4f   (rel err %.2e)\n', ...
        CoC_K0, C_T2_ATM, abs(CoC_K0 - C_T2_ATM)/C_T2_ATM);

%% ====================== 4a Put-on-Put (PoP) ===========================
rng(seed);
[PoP_price, PoP_IC] = price_POP_GL_MC(T1, T2, alpha_GL, beta_GL, ...
                                      sigma_T1, sigma_T2, Nsim, ...
                                      forward(iT2), discount_factor(iT1), ...
                                      discount_factor(iT2), N_grid, K1);

% --- Sanity check (calcolo SEPARATO, K1 "deep" grande) ----------------
% Per il PoP il collasso utile NON e' K1=0 (darebbe 0), ma K1 grande:
% se K1 supera sempre P(T1), allora max(K1 - P, 0) = K1 - P e il max sparisce.
%   PoP_deep = B(0,T1)*K1 - B(0,T1)*E[P(T1)]
% e per tower property + martingalita' B(0,T1)*E[P(T1)] e' la vanilla put
% ATM-forward a T2:  put_T2_ATM = B(0,T2)*sigma_t(T2)*sqrt(T2)*G0
% (centratura => E[z^-]=E[z^+]=G0). 
% Quindi:PoP_deep_ref = B(0,T1)*K1_deep - put_T2_ATM
% Scegliamo K1_deep molto sopra il valore atteso della put interna a T1
% (la put interna a T1 vale in media ~ put ATM a T2 / B(0,T1)).
put_T2_ATM = discount_factor(iT2) * sigma_t(iT2) * sqrt(T2) * G0;
K1_deep    = 50 * put_T2_ATM / discount_factor(iT1);   % abbondante margine

rng(seed);
[PoP_deep, ~] = price_POP_GL_MC(T1, T2, alpha_GL, beta_GL, ...
                                sigma_T1, sigma_T2, Nsim, ...
                                forward(iT2), discount_factor(iT1), ...
                                discount_factor(iT2), N_grid, K1_deep);
PoP_deep_ref = discount_factor(iT1) * K1_deep - put_T2_ATM;

fprintf('\n--- 4a PoP GL (K1=%.2f, K2=F(t0,T2)=%.4f) ---\n', K1, forward(iT2));
fprintf('MC price       : %.4f   95%% CI [%.4f, %.4f]\n', PoP_price, PoP_IC(1), PoP_IC(2));
fprintf('Sanity (deep)  : MC %.4f  vs  ref %.4f   (rel err %.2e)\n', ...
        PoP_deep, PoP_deep_ref, abs(PoP_deep - PoP_deep_ref)/abs(PoP_deep_ref));

%% ====================== 4b Chooser ====================================
rng(seed);
[ch_price, ch_IC] = price_chooser_GL_MC(T1, T2, alpha_GL, beta_GL, ...
                                        sigma_T1, sigma_T2, Nsim, ...
                                        forward(iT2), discount_factor(iT2), N_grid);

% --- Sanity check analitico (Stochastic-Reset) ------------------------
% Chooser ATM = B(0,T2)*G0*[sigma_t(T2)*sqrt(T2) + sigma_t(T1)*sqrt(T1)].
% Derivazione: max(C,P) = C + B(T1,T2)*(K - F(T1,T2))^+ ; (I) -> call ATM
% a T2, (II) -> put ATM su f_{T1}. Centratura (media 0) => E[z^-]=E[z^+]=G0.
% USA sigma_t per coerenza con i pricer.
ch_analytic = discount_factor(iT2) * G0 * ...
              ( sigma_t(iT2)*sqrt(T2) + sigma_t(iT1)*sqrt(T1) );

fprintf('\n--- 4b Chooser GL (K=F(t0,T2)=%.4f, T1=%.2fy, T2=%.2fy) ---\n', ...
        forward(iT2), T1, T2);
fprintf('MC price       : %.4f   95%% CI [%.4f, %.4f]\n', ch_price, ch_IC(1), ch_IC(2));
fprintf('Analytic (SR)  : %.4f                            (rel err %.2e)\n', ...
        ch_analytic, abs(ch_price - ch_analytic)/ch_analytic);

%% ====================== Riepilogo =====================================
fprintf('\n========== RIEPILOGO PREZZI (K1=%.2f) ==========\n', K1);
fprintf('CoC      : %.4f\n', CoC_price);
fprintf('PoP      : %.4f\n', PoP_price);
fprintf('Chooser  : %.4f\n', ch_price);