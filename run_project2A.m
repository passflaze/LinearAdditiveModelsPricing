% data loading
clear,clc
addpath("Utilities/");
addpath("ex2/");
addpath("ex3_utilities/");
addpath("ex4/");

callpath="Data/datacalls";
putpath="Data/dataputs";
expiryFile = "Data/Expiries_Futures.txt";
valueDate = datetime(2020,06,02);

[strikes, calls, puts, expiries] = readData(callpath, putpath, valueDate, expiryFile);

%% bootstrap: synthetic discount factor and forward, one per maturity
nT = numel(expiries);
discount_factor = zeros(nT,1);
forward = zeros(nT,1);
R2 = zeros(nT,1);
call_atm = zeros(nT,1);

for k = 1:nT
    [discount_factor(k), forward(k), R2(k)] = bootstrap(puts(k,:), calls(k,:), strikes);
    call_atm(k) = callATM(calls(k,:), puts(k,:), strikes, forward(k), discount_factor(k));
end

act_365 = 3; % actual/365 day count convention for time to maturity
yf = yearfrac(valueDate, expiries, act_365);
sigma_atm = sigmaATM(call_atm, discount_factor, yf, expiries);

%% OTM filter -> pooled calibration of (eta, k)
% Moneyness range [-30$, 30$] as in Baviera & Massaria (2026), Table 3.
% Following the paper (Eq. 20) the SSE is on full $ prices: keep cMkt in $ and
% pass per-point norm = B * sigma_ATM * sqrt(t) so the objective can rescale
% the normalized model output to dollars.
xMax = 30;

chi_cell      = cell(nT, 1);
cMkt_cell     = cell(nT, 1);
normFact_cell = cell(nT, 1);

for iT = 1:nT
    F_k  = forward(iT);
    B_k  = discount_factor(iT);
    t_k  = yf(iT);
    s_k  = sigma_atm(iT);
    norm_iT = B_k * s_k * sqrt(t_k);

    % OTM calls: 0 < K - F <= xMax, quoted only (drop untraded NaN prices)
    otm_c = (strikes > F_k) & (strikes <= F_k + xMax) & isfinite(calls(iT,:));
    K_c   = strikes(otm_c)';
    c_c   = calls(iT, otm_c)';

    % OTM puts: -xMax <= K - F < 0 -> converted to call via put-call parity
    % C = P + B*(F-K)
    otm_p = (strikes < F_k) & (strikes >= F_k - xMax) & isfinite(puts(iT,:));
    K_p   = strikes(otm_p)';
    c_p   = puts(iT, otm_p)' + B_k * (F_k - K_p);

    nPts = numel(K_p) + numel(K_c);
    chi_cell{iT}      = ([ K_p; K_c ] - F_k) / (s_k * sqrt(t_k));
    cMkt_cell{iT}     = [ c_p; c_c ];
    normFact_cell{iT} = norm_iT * ones(nPts, 1);
end

chi_all      = vertcat(chi_cell{:});
cMkt_all     = vertcat(cMkt_cell{:});
normFact_all = vertcat(normFact_cell{:});

%% calibrate Additive Bachelier (global: eta and k constant across maturities)
% paper Eq. (20): SSE on dollar prices.
[kAB, eta, sigma_t, rMSE] = calibrateAB(chi_all, cMkt_all, normFact_all, sigma_atm);

%% verify the fit visually in normalized space (cleaner across maturities)
% G_model = G(chi; eta, k) is the normalized model price (call_AB_FFT).
% market normalized price: cNorm = cMkt / (B * sigma_ATM * sqrt(t))
cNorm_all = cMkt_all ./ normFact_all;

chi_grid = linspace(min(chi_all), max(chi_all), 200)';
G_model  = call_AB_FFT(chi_grid, kAB, eta);

figure; hold on
scatter(chi_all, cNorm_all, 12, 'filled')        % market
plot(chi_grid, G_model, 'r-', 'LineWidth', 1.5)  % model
xlabel('\chi'); ylabel('$\mathcal{G}(\chi)$','Interpreter','latex')
legend('market','AB model'); grid on
title(sprintf('AB calibration (paper Eq. 20): k = %.4f, \\eta = %.4f, rMSE = %.4f $', ...
              kAB, eta, rMSE));

%% 3b-3c: conditional CDF and MC simulation of the log-price increment T1 -> T2
% T1 = 6m, T2 = 1y: closest available maturities idx 3 (5.52m) and idx 5 (11.47m).
seed = 1234;
rng(seed);  % for reproducibility
iT1 = 3;  iT2 = 5;
T1 = yf(iT1);  T2 = yf(iT2);
sigma_T1 = sigma_t(iT1);  sigma_T2 = sigma_t(iT2);

% 3b: conditional CDF via Lewis-FFT (Baviera-Manzoni 2026, eq. 13-15)
x_grid = linspace(-40,40,300)';


cdf    = ccdf_AB_FFT(eta, kAB, T1, T2, sigma_T1, sigma_T2, x_grid, 0);

% 3c: simulate increments via inverse-CDF + exponential tail extrapolation
Nsim = 1e5;
Z    = sample_from_cdf(x_grid, cdf, Nsim);

% visual check: MC histogram/CDF vs FFT theoretical
plot_mc_check(Z, x_grid, cdf, T1, T2);

% 3d: forward-start option pricing  K=1 -> payoff = max(Z_{2|1}, 0)
[price, IC] = price_AB_MC(T1, T2, kAB, eta, sigma_T1, sigma_T2, Nsim, ...
                           forward(iT2), discount_factor(iT2), x_grid, 1);

% FFT reference: E[max(Z_{2|1},0)] = integral_0^inf (1-F_cond(x)) dx
% x_grid may not contain x=0 exactly, so we interpolate F(0) to avoid
% truncating the integral over [0, first positive grid point].
F0       = interp1(x_grid, cdf, 0, 'spline');
mask_pos = x_grid > 0;
price_fft = discount_factor(iT2) * trapz([0; x_grid(mask_pos)], ...
                                          1 - [F0; cdf(mask_pos)]);

fprintf('--- Forward-start option (K=1, T1=%.2fy, T2=%.2fy) ---\n', T1, T2);
fprintf('FFT reference : %.4f\n', price_fft);
fprintf('MC estimate   : %.4f  95%% CI [%.4f, %.4f]\n', price, IC(1), IC(2));

%% 4a Call-on-Call (CoC) via Lewis-FFT + MC
% Inner strike K2 = F(t0,T2); CoC payoff at T1 = max( C(T1; K2, T2) - K1, 0 ).
K1 = 1;   % outer strike (placeholder, set per project spec)

var_T1   = sigma_T1^2 * T1 * (1 + eta^2 * kAB);
var_T2   = sigma_T2^2 * T2 * (1 + eta^2 * kAB);
x_grid   = linspace(-8*sqrt(var_T2 - var_T1), 8*sqrt(var_T2 - var_T1), 300)';

[CoC_price, CoC_IC] = price_COC_AB_MC(T1, T2, kAB, eta, sigma_T1, sigma_T2, Nsim, ...
                                      forward(iT2), discount_factor(iT1), ...
                                      discount_factor(iT2), x_grid, K1);

% Sanity: K1=0 collapses the outer max -> vanilla ATM-forward Call at T2.
[CoC_K0, ~] = price_COC_AB_MC(T1, T2, kAB, eta, sigma_T1, sigma_T2, Nsim, ...
                              forward(iT2), discount_factor(iT1), ...
                              discount_factor(iT2), x_grid, 0);
G0       = call_AB_FFT(0, kAB, eta);                          % = 1/sqrt(2*pi)
C_T2_ATM = discount_factor(iT2) * sigma_atm(iT2) * sqrt(T2) * G0;

fprintf('\n--- 4a CoC (K1=%.2f, K2=F(t0,T2)=%.4f, T1=%.2fy, T2=%.2fy) ---\n', ...
        K1, forward(iT2), T1, T2);
fprintf('MC price       : %.4f   95%% CI [%.4f, %.4f]\n', CoC_price, CoC_IC(1), CoC_IC(2));
fprintf('Sanity (K1=0)  : MC %.4f  vs  FFT Call_T2_ATM %.4f   (rel err %.2e)\n', ...
        CoC_K0, C_T2_ATM, abs(CoC_K0 - C_T2_ATM)/C_T2_ATM);

%% 4b Chooser via Lewis-FFT + put-call parity + MC
% Choice at T1 between Call/Put with K = F(t0,T2), maturity T2.
Nsim =1e6;
[ch_price, ch_IC] = price_chooser_AB_MC(T1, T2, kAB, eta, sigma_T1, sigma_T2, Nsim, ...
                                        forward(iT2), discount_factor(iT2), x_grid);

% Sanity: max(C,P) = C + B(T1,T2)*(K-F(T1,T2))+ and martingality of f_t give
%   Chooser = B(0,T2) * G(0) * [sigma_ATM(T2)*sqrt(T2) + sigma_ATM(T1)*sqrt(T1)]
% (smile term drops at ATM-forward both at T1 and T2).
ch_analytic = discount_factor(iT2) * G0 * ...
              ( sigma_atm(iT2)*sqrt(T2) + sigma_atm(iT1)*sqrt(T1) );

fprintf('\n--- 4b Chooser (K=F(t0,T2)=%.4f, T1=%.2fy, T2=%.2fy) ---\n', ...
        forward(iT2), T1, T2);
fprintf('MC price       : %.4f   95%% CI [%.4f, %.4f]\n', ch_price, ch_IC(1), ch_IC(2));
fprintf('Analytic (SR)  : %.4f                            (rel err %.2e)\n', ...
        ch_analytic, abs(ch_price - ch_analytic)/ch_analytic);
