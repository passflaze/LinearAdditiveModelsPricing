% data loading
clear,clc
addpath("Utilities/");
addpath("ex2/");

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

%% OTM filter + normalization → global calibration of eta and k
% moneyness range [-30$, 30$] as in Baviera & Massaria (2026), Table 3
xMax = 30;

chi_cell   = cell(nT, 1);
cNorm_cell = cell(nT, 1);

for iT = 1:nT
    F_k  = forward(iT);
    B_k  = discount_factor(iT);
    t_k  = yf(iT);
    s_k  = sigma_atm(iT);
    norm = B_k * s_k * sqrt(t_k);   % normalization factor

    % OTM calls: 0 < K - F <= xMax, quoted only (drop untraded NaN prices)
    otm_c = (strikes > F_k) & (strikes <= F_k + xMax) & isfinite(calls(iT,:));
    K_c   = strikes(otm_c)';
    c_c   = calls(iT, otm_c)' / norm;

    % OTM puts: -xMax <= K - F < 0 → converted to calls via put-call parity
    % C = P + B*(F-K)  =>  C_norm = (P + B*(F-K)) / norm
    otm_p = (strikes < F_k) & (strikes >= F_k - xMax) & isfinite(puts(iT,:));
    K_p   = strikes(otm_p)';
    c_p   = ( puts(iT, otm_p)' + B_k * (F_k - K_p) ) / norm;

    chi_cell{iT}   = ([ K_p; K_c ] - F_k) / (s_k * sqrt(t_k));
    cNorm_cell{iT} = [ c_p; c_c ];
end

chi_all   = vertcat(chi_cell{:});
cNorm_all = vertcat(cNorm_cell{:});

%% calibrate Additive Bachelier (global: eta and k constant across maturities)
[kAB, eta, sigma_t, MSE] = calibrateAB(chi_all, cNorm_all, sigma_atm)

%% verify the fit visually: scatter of market vs model price on the chi grid
chi_grid = linspace(min(chi_all), max(chi_all), 200)';
G_model  = call_AB_FFT(chi_grid, kAB, eta);

figure; hold on
scatter(chi_all, cNorm_all, 12, 'filled')      % mercato
plot(chi_grid, G_model, 'r-', 'LineWidth', 1.5) % modello
xlabel('\chi'); ylabel('$\hat{G}(\chi)$','Interpreter','latex')
legend('mercato','modello AB'); grid on


