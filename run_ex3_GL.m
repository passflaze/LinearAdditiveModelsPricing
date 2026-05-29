% data loading
clear, clc
addpath("Utilities/");
addpath("Functions/");
addpath("ex3_GL/");

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

%% RUN EX 3 - GL model
% CF dell'incremento via additivita'

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
s = yearfrac(valueDate, valueDate + calmonths(6), 3);
t = yearfrac(valueDate, valueDate + calyears(1),  3);
sigma_s = sigma_ATM(2) / I0_opt;
sigma_t = sigma_ATM(4) / I0_opt;


[CDF_grid, zk] = lewis_fft_cdf_old(@cf_increment_GL, alpha_GL, beta_GL, sigma_s, s, sigma_t, t);

fprintf('CDF estremi: P(x_1)=%.4e   1-P(x_N)=%.4e\n', ...
        CDF_grid(1), 1 - CDF_grid(end));

%% Plot CDF
figure('Name', 'Check CDF FFT - Metodo Splicing', 'Position', [200, 200, 700, 500]);
plot(zk, CDF_grid, 'b-', 'LineWidth', 2);
hold on;
yline(0, 'k--', 'Alpha', 0.3);
yline(1, 'k--', 'Alpha', 0.3);
xline(0, 'k--', 'Alpha', 0.3);
title('CDF Incremento - Splicing Doppio Shift');
xlabel('z');
ylabel('Cumulative Probability');
grid on;
hold off;

fprintf('CDF estremi: P(z_1)=%.4e   1-P(z_N)=%.4e\n', ...
        CDF_grid(1), 1 - CDF_grid(end));

%% CHECK: VALIDAZIONE ALGORITMO FFT SU NORMALE STANDARD
fprintf('\n--- Validazione Algoritmo FFT (Gaussiana) ---\n');

% CF della normale standard con firma compatibile (ignora i parametri extra)
cf_gauss = @(u, alpha, beta, sigma_s, s, sigma_t, t) exp(-0.5 .* u.^2);

% Parametri scelti per forzare a_neg = -1, a_pos = +1
% (a_neg = -alpha/(2*sigma_t*sqrt(t)),  a_pos = beta/(2*sigma_t*sqrt(t)))
alpha_g   = 2;
beta_g    = 2;
sigma_s_g = 1;   % irrilevante (ignorato dalla cf_gauss)
s_g       = 1;   % irrilevante
sigma_t_g = 1;
t_g       = 1;

% Chiamata alla funzione
[CDF_gauss, zk_g] = lewis_fft_cdf_old(cf_gauss, alpha_g, beta_g, ...
                                  sigma_s_g, s_g, sigma_t_g, t_g);

% CDF analitica
CDF_analytical = normcdf(zk_g);

% Confronto numerico
max_err = max(abs(CDF_gauss - CDF_analytical));
fprintf('Errore Massimo (FFT vs Analitica Normale): %.2e\n', max_err);

% Plot
figure('Name', 'Check FFT su Gaussiana', 'Position', [250, 250, 700, 500]);
plot(zk_g, CDF_gauss,      'b-',  'LineWidth', 3, 'DisplayName', 'FFT Inversion');
hold on;
plot(zk_g, CDF_analytical, 'r--', 'LineWidth', 2, 'DisplayName', 'Analytical normcdf');
title('Validazione Algoritmo: FFT vs Analytical Gaussian CDF');
xlabel('z');
ylabel('Cumulative Probability');
legend('Location', 'northwest');
xlim([-5, 5]);
grid on;
hold off;

%% SIMULAZIONE

[X_st, x_grid, CDF_clean] = simulate_increments(zk, CDF_grid, 1e5, 42);

%% PLOT FINALI CHECK
figure('Name', 'Lewis-FFT-S', 'Position', [200 200 1200 400]);

subplot(1,3,1);
plot(zk, CDF_grid, 'b-', 'LineWidth', 1.5); hold on;
yline(0, 'k--', 'Alpha', 0.3);
yline(1, 'k--', 'Alpha', 0.3);
xline(x_grid(1),   'r--');
xline(x_grid(end), 'r--');
xlim([x_grid(1) - 5, x_grid(end) + 5]);
xlabel('x');
ylabel('$\hat{P}(x)$', 'Interpreter', 'latex');
title('CDF spliced (zona utile in rosso)');
grid on; hold off;

subplot(1,3,2);
[f_emp, x_emp] = ecdf(X_st);
plot(x_emp, f_emp, 'r-', 'LineWidth', 1.5); hold on;
plot(x_grid, CDF_clean, 'b--', 'LineWidth', 1);
xlim([x_grid(1), x_grid(end)]);
xlabel('x'); ylabel('P(x)');
legend('empirica', 'target', 'Location', 'best');
title('CDF empirica vs target');
grid on; hold off;

subplot(1,3,3);
histogram(X_st, 100, 'Normalization', 'pdf', ...
          'FaceColor', '#4DBEEE', 'EdgeColor', 'none');
xlim([x_grid(1), x_grid(end)]);
xlabel('x'); ylabel('pdf');
title('PDF empirica simulata');
grid on;

%% prezzo MC vs prezzo ANALITICO 

strike =1;
price_MC = price_forward_start_MC(alpha_GL, beta_GL, sigma_s, s, sigma_t, t, ...
                            forward(4), discount_factor(4), strike, 1e5, 42);

%prezzo analitico
cdf_fun = @(x) interp1(x_grid, CDF_clean, x, 'spline', 'extrap');
integrand = @(x) 1 - cdf_fun(x);
expectation_analytical = integral(integrand, 0, x_grid(end));
price_analytical = discount_factor(4) * expectation_analytical;

absolute_error = abs(price_MC - price_analytical);
percentage_error = (absolute_error / price_analytical) * 100;

fprintf('FFT reference (Analitico) : %.4f\n', price_analytical);
fprintf('MC estimate (1e6 paths)   : %.4f\n', price_MC);
fprintf('Errore Assoluto                  : %.2e\n', absolute_error);
fprintf('Errore Percentuale               : %.4f %%\n', percentage_error);


