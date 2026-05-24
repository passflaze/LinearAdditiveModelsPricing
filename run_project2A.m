% data loading
clear,clc
addpath("Utilities/");
addpath("ex2/");
addpath("ex3_utilities/");

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

%% 3b: conditional CDF of the log-price increment between T1 and T2 via Lewis-FFT
% Point 3c requires T1 = 6m, T2 = 1y -> closest available maturities are
% idx 3 (DEC20, 5.52m) and idx 5 (JUN21, 11.47m).
iT1 = 3;
iT2 = 5;
T1 = yf(iT1);
T2 = yf(iT2);
sigma_T1 = sigma_t(iT1);
sigma_T2 = sigma_t(iT2);

x_grid = linspace(-40, 40, 300)';        
cdf = ccdf_AB_FFT(eta, kAB, T1, T2, sigma_T1, sigma_T2, x_grid);

idx_b = find(cdf >= 0, 1, 'first'); % first idx where cdf >= 0
idx_e = find(cdf <= 1, 1, 'last'); % last idx where cdf <= 1 
x_valid = x_grid(idx_b:idx_e);
cdf_valid = cdf(idx_b:idx_e);
fprintf('Tail probability mass = %.2e : less than 1e-4', max(cdf_valid(1), 1 - cdf_valid(end)) );

Nsim = 1e5;
U    = rand(Nsim, 1);
Z    = interp1(cdf_valid, x_valid, U, 'spline');

xb = x_valid(1);  xe = x_valid(end);
Pb = cdf_valid(1); Pe = cdf_valid(end);
lam_m = (log(cdf_valid(2)) - log(Pb)) / (x_valid(2) - xb);
lam_p = (log(1 - cdf_valid(end-1)) - log(1 - Pe)) / (x_valid(end-1) - xe);

% per coda sinistra (U < Pb):  U = Pb * exp(lam_m * (z - xb))  =>  z = xb + log(U/Pb)/lam_m
% per coda destra (U > Pe):    1-U = (1-Pe) * exp(-lam_p * (z - xe)) => z = xe - log((1-U)/(1-Pe))/lam_p
left  = U < Pb;
right = U > Pe;
Z(left)  = xb + log(U(left)/Pb) / lam_m;
Z(right) = xe - log((1 - U(right))/(1 - Pe)) / lam_p;


%% visualizzazione incrementi simulati T1 -> T2
figure;

% pannello 1: istogramma simulato vs PDF teorica (da CDF per diff finite)
subplot(1,2,1)
histogram(Z, 80, 'Normalization', 'pdf', 'FaceAlpha', 0.5)
hold on
pdf_approx = diff(cdf_valid) ./ diff(x_valid);   % f ≈ dF/dx
x_mid      = (x_valid(1:end-1) + x_valid(2:end)) / 2;
plot(x_mid, pdf_approx, 'r-', 'LineWidth', 1.5)
xlabel('x = X_{T_2} - X_{T_1} [$]')
ylabel('pdf')
legend('MC', 'teorica (FFT)')
title(sprintf('Incremento AB: T_1=%.2fy \\to T_2=%.2fy', T1, T2))
grid on

% pannello 2: CDF empirica vs CDF teorica
subplot(1,2,2)
[f_emp, x_emp] = ecdf(Z);
plot(x_emp, f_emp, 'b-', 'LineWidth', 1.2); hold on
plot(x_valid, cdf_valid, 'r--', 'LineWidth', 1.5)
xlabel('x [$]'); ylabel('F(x)')
legend('MC empirica', 'FFT teorica')
xlim([-30 30]); grid on
title('CDF: MC vs modello')














