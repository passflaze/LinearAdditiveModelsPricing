function [params, market] = run_ex2(opts)
% RUN_EX2  Calibrate the AB, MA, and GL models to the WTI vol surface.
%   Performs curve bootstrapping, ATM vol extraction, OTM surface filtering,
%   and model calibration. Returns reusable structs consumed by run_ex3/4/6.
%
% INPUTS:
%   opts - (struct, optional) all fields have defaults:
%     .callpath    - calls folder       (default "Data/datacalls")
%     .putpath     - puts folder        (default "Data/dataputs")
%     .expiryFile  - expiry file path   (default "Data/Expiries_Futures.txt")
%     .valueDate   - value date         (default 2020-06-02)
%     .x_min       - dollar-moneyness lower bound  (default -30)
%     .x_max       - dollar-moneyness upper bound  (default  30)
%     .M           - FFT grid exponent for GL      (default 15)
%     .dz          - FFT grid step for GL          (default 2.5e-3)
%     .verbose     - print full report             (default false)
%     .plot        - draw term-structure figures   (default false)
%
% OUTPUTS:
%   params - (struct) calibrated parameters (column vectors):
%              .AB -> [k; eta]  .MA -> [alpha; beta]  .GL -> [alpha; beta]
%   market - (struct) market and support data (column vectors):
%              .strikes .calls .puts .expiries .valueDate
%              .discount_factor .forward .R2
%              .c_ATM .sigma_ATM .yf .moneyness_modified .c_mkt_calibration

addpath("Calibration/Calibration_MA/");
addpath("Calibration/Calibration_AB/");
addpath("Calibration/Calibration_GL/");

if nargin < 1 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'callpath'),   opts.callpath   = "Data/datacalls";            end
if ~isfield(opts, 'putpath'),    opts.putpath    = "Data/dataputs";             end
if ~isfield(opts, 'expiryFile'), opts.expiryFile = "Data/Expiries_Futures.txt"; end
if ~isfield(opts, 'valueDate'),  opts.valueDate  = datetime(2020, 06, 02);      end
if ~isfield(opts, 'x_min'),      opts.x_min      = -30;                         end
if ~isfield(opts, 'x_max'),      opts.x_max      =  30;                         end
if ~isfield(opts, 'M'),          opts.M          = 15;                          end
if ~isfield(opts, 'dz'),         opts.dz         = 2.5e-3;                       end
if ~isfield(opts, 'verbose'),    opts.verbose    = false;                        end
if ~isfield(opts, 'plot'),       opts.plot       = false;                       end

vb = @(varargin) verbose_print(opts.verbose, varargin{:});
fprintf('=========================================================================\n');
fprintf('             VOLATILITY SURFACE CALIBRATION ENGINE                       \n');
fprintf('=========================================================================\n\n');

%% --- Data loading and bootstrapping ---
vb('STEP 1: Loading data and bootstrapping curve...\n');

[strikes, calls, puts, expiries] = readData(opts.callpath, opts.putpath, ...
                                             opts.valueDate, opts.expiryFile);

nT = numel(expiries);
discount_factor = zeros(nT, 1);
forward         = zeros(nT, 1);
R2              = zeros(nT, 1);

for k = 1:nT
    [discount_factor(k), forward(k), R2(k)] = bootstrap(puts(k,:), calls(k,:), strikes);
end

vb("\n--- Bootstrap Results (Value Date: %s) ---\n", string(opts.valueDate, "yyyy-MM-dd"));
vb("%-12s  %10s  %10s  %8s\n", "Expiry", "D(t,T)", "F(t,T)", "R^2");
for k = 1:nT
    vb("%-12s  %10.6f  %10.4f  %8.4f\n", ...
        string(expiries(k), "yyyy-MM-dd"), ...
        discount_factor(k), forward(k), R2(k));
end
vb('-------------------------------------------------------------------------\n\n');

yf = yearfrac(opts.valueDate, expiries, 3);
zero_rates = -log(discount_factor) ./ yf;

%% --- ATM vol calibration and moneyness generation ---
vb('STEP 2: Calibrating ATM Volatility and filtering surface...\n');
c_ATM = zeros(length(forward), 1);
for i = 1:length(forward)
    current_calls = calls(i, :);
    current_puts  = puts(i, :);
    c_ATM(i) = callATM(current_calls, current_puts, strikes, forward(i), discount_factor(i));
end

sigma_ATM = sigmaATM(c_ATM, discount_factor, yf, expiries, opts.plot);
if opts.verbose
    check_term_structure(sigma_ATM, yf, expiries);
else
    evalc('check_term_structure(sigma_ATM, yf, expiries);');   % silence diagnostic
end

x_min = opts.x_min;
x_max = opts.x_max;
[moneyness_modified, c_mkt_calibration] = moneyness_generator( ...
    forward, strikes, calls, puts, sigma_ATM, yf, discount_factor, ...
    x_min, x_max);
vb('  -> Surface filtered on dollar moneyness x in [%g, %g] $ and prepared for optimization.\n\n', ...
    x_min, x_max);

%% --- AB model calibration ---
vb('STEP 3: Calibrating Additive Bachelier (AB) Model via fmincon...\n');

if opts.verbose, ab_disp = 'iter'; else, ab_disp = 'off'; end
options_AB = optimoptions('fmincon', ...
    'Display', ab_disp, ...
    'Algorithm', 'interior-point', ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance',       1e-10, ...
    'MaxFunctionEvaluations', 5000);

x0_AB = [1.0, 0.2];
lb_AB = [1e-3, -1.5];
ub_AB = [5.0,   1.5];
obj_fun_AB = @(x) objective_function_AB(x, discount_factor, yf, sigma_ATM, ...
                       moneyness_modified, c_mkt_calibration);
                   
[x_opt_AB, fval_AB, exitflag_AB] = ...
    fmincon(obj_fun_AB, x0_AB, [], [], [], [], lb_AB, ub_AB, [], options_AB);
eta_AB = x_opt_AB(2);
k_AB   = x_opt_AB(1);

vb('\n  -> Optimal parameters found: k = %.6f, eta = %.6f.\n', k_AB, eta_AB);
vb('  -> AB Calibration completed (exitflag = %d, SSE = %.6g).\n\n', exitflag_AB, fval_AB);

%% --- MA model calibration ---
vb('STEP 4: Calibrating Minimal Additive (MA) Model via fminbnd...\n');

ALPHA_FIX = 1.0;
lb_beta = 0.05;
ub_beta = 40;

if opts.verbose, ma_disp = 'iter'; else, ma_disp = 'off'; end
options = optimoptions('fmincon', ...
    'Display', ma_disp, ...
    'Algorithm', 'interior-point', ...
    'FunctionTolerance',   1e-12, ...
    'StepTolerance',       1e-12, ...
    'OptimalityTolerance', 1e-10, ...
    'MaxFunctionEvaluations', 2000);

obj_fun_MA_1d = @(b) objective_function_MA([ALPHA_FIX, b], ...
                       discount_factor, yf, sigma_ATM, ...
                       moneyness_modified, c_mkt_calibration);

opts_MA_1d = optimset('TolX', 1e-10, 'Display', ma_disp);
[beta_opt, fval_MA, exitflag_MA] = fminbnd(obj_fun_MA_1d, lb_beta, ub_beta, opts_MA_1d);

alpha_MA = ALPHA_FIX;
beta_MA  = beta_opt;

vb('\n  -> MA Calibration completed (exitflag = %d, SSE = %.6g).\n', exitflag_MA, fval_MA);
vb('  -> Fixed scale: alpha = %.6f (gauge); calibrated beta = %.6f.\n', alpha_MA, beta_MA);
vb('  -> Asymmetry ratio beta/alpha = %.6f.\n\n', beta_MA/alpha_MA);

% Heavy MA skew diagnostic: prints a report AND draws a figure. Run it only
% when the full report is requested (verbose) or plots are on; keep the console
% clean otherwise (and skip the stray figure when plot=false).
if opts.verbose
    check_skew_MA(alpha_MA, beta_MA, discount_factor, yf, ...
        sigma_ATM, moneyness_modified, c_mkt_calibration, expiries);
elseif opts.plot
    evalc(['check_skew_MA(alpha_MA, beta_MA, discount_factor, yf, ', ...
           'sigma_ATM, moneyness_modified, c_mkt_calibration, expiries);']);
end

%% --- GL model calibration ---
vb('STEP 5: Calibrating Generalized Laplace (GL) Model via fmincon...\n');

M  = opts.M;
dz = opts.dz;
nonlcon_TS = @(x) term_structure_nonlcon(x, sigma_ATM, yf);

x0_GL = [0.82, 1.22];
lb_GL = [0.05, 0.05];
ub_GL = [50,  50];
obj_fun_GL = @(x) objective_function_GL(x, discount_factor, yf, sigma_ATM, moneyness_modified, c_mkt_calibration, M, dz);

[x_opt_GL, fval_GL, exitflag_GL] = fmincon(obj_fun_GL, x0_GL, [], [], [], [], lb_GL, ub_GL, nonlcon_TS, options);
alpha_GL = x_opt_GL(1);
beta_GL  = x_opt_GL(2);

vb('\n  -> Optimal parameters found: alpha = %.6f, beta = %.6f.\n', alpha_GL, beta_GL);
vb('  -> GL Calibration completed (exitflag = %d, SSE = %.6g).\n\n', exitflag_GL, fval_GL);

%% --- Diagnostics and output packing ---
if opts.verbose
    print_diagnostics(k_AB, eta_AB, fval_AB, alpha_MA, beta_MA, fval_MA, ...
        alpha_GL, beta_GL, fval_GL, M, dz, discount_factor, yf, sigma_ATM, ...
        moneyness_modified, c_mkt_calibration, expiries, nT, opts.plot);
else
    % Compact report (verbose=false): discount factors + calibrated parameters.
    fprintf('\n--- Discount Factors (Value Date: %s) ---\n', string(opts.valueDate, "yyyy-MM-dd"));
    fprintf('%-12s  %12s  %12s\n', 'Expiry', 'D(t,T)', 'F(t,T)');
    for kk = 1:nT
        fprintf('%-12s  %12.6f  %12.4f\n', string(expiries(kk), "yyyy-MM-dd"), ...
                discount_factor(kk), forward(kk));
    end
    fprintf('\n--- Calibrated Parameters ---\n');
    fprintf('  AB : k     = %9.6f | eta  = %9.6f   (SSE = %.4g)\n', k_AB, eta_AB, fval_AB);
    fprintf('  MA : alpha = %9.6f | beta = %9.6f   (SSE = %.4g)\n', alpha_MA, beta_MA, fval_MA);
    fprintf('  GL : alpha = %9.6f | beta = %9.6f   (SSE = %.4g)\n\n', alpha_GL, beta_GL, fval_GL);
end

params = struct();
params.AB = [k_AB;     eta_AB ];
params.MA = [alpha_MA; beta_MA];
params.GL = [alpha_GL; beta_GL];

if nargout > 1
    market = struct();
    market.strikes            = strikes;
    market.calls              = calls;
    market.puts               = puts;
    market.expiries           = expiries;
    market.valueDate          = opts.valueDate;
    market.discount_factor    = discount_factor;
    market.forward            = forward;
    market.R2                 = R2;
    market.c_ATM              = c_ATM;
    market.sigma_ATM          = sigma_ATM;
    market.yf                 = yf;
    market.moneyness_modified = moneyness_modified;
    market.c_mkt_calibration  = c_mkt_calibration;
    market.zero_rates         = zero_rates;
end

%% --- Plotting suite ---
if opts.plot
    % --- 7.1: Discount Factors and Zero Rates ---
    figure('Name', 'Discount Factors and Zero Rates', 'Color', 'w');
    yyaxis left
    plot(expiries, discount_factor, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    ylabel('Discount Factor');
    ax = gca;
    ax.YColor = 'r'; 
    yyaxis right
    plot(expiries, zero_rates, 'b-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    ylabel('Zero Rate (Continuously Compounded)');
    ax.YColor = 'b'; 
    xlabel('Expiry Date');
    title('Market Discount Factors and Zero Rates');
    grid on;
    if isdatetime(expiries)
        xtickformat('yyyy-MM-dd');
    end

    % --- 7.2: Implied-volatility term structure at a fixed wing ---
    x0 = 20;                                  % fixed dollar moneyness K - F ($), OTM call wing
    chi_w = x0 ./ (sigma_ATM .* sqrt(yf));    % (M x 1) per-maturity normalized moneyness at x0

    % model call prices at (chi_w, each maturity) -> (M x 1) each
    c_AB_w = price_AB(params.AB, discount_factor, yf, sigma_ATM, chi_w);
    c_MA_w = price_MA(params.MA, discount_factor, yf, sigma_ATM, chi_w);
    c_GL_w = price_GL(params.GL(1), params.GL(2), M, dz, ...
                      discount_factor, sigma_ATM, yf, chi_w);

    % invert to Bachelier implied vol (F - K = -x0 for every maturity)
    fmk_w     = -x0 * ones(size(yf));
    sigIV_AB  = bachelier_iv(c_AB_w, discount_factor, fmk_w, yf);
    sigIV_MA  = bachelier_iv(c_MA_w, discount_factor, fmk_w, yf);
    sigIV_GL  = bachelier_iv(c_GL_w, discount_factor, fmk_w, yf);

    figure('Name', 'Implied-vol term structure at a fixed wing', 'Color', 'w');
    hold on;
    plot(expiries, sigma_ATM, 'k--', 'LineWidth', 1.5, 'DisplayName', '\sigma_{ATM} (x=0)');
    plot(expiries, sigIV_AB,  'b-',  'LineWidth', 1.5, 'DisplayName', sprintf('AB  (x=%g$)', x0));
    plot(expiries, sigIV_MA,  'r-',  'LineWidth', 1.5, 'DisplayName', sprintf('MA  (x=%g$)', x0));
    plot(expiries, sigIV_GL,  'g-',  'LineWidth', 1.5, 'DisplayName', sprintf('GL  (x=%g$)', x0));
    xlabel('Maturity');
    ylabel('Bachelier implied volatility ($)');
    title(sprintf('Implied-vol term structure: ATM vs models at x = K-F = %g$', x0));
    legend('Location', 'best');
    grid on;
    if isdatetime(expiries), xtickformat('yyyy-MM'); end
    hold off;

    % --- 7.3: Distributions and Implied Volatility Smile (Bachelier) ---
    plot_distributions(eta_AB, k_AB, alpha_MA, beta_MA, alpha_GL, beta_GL);
    
    plot_iv_AB(k_AB, eta_AB, discount_factor, yf, sigma_ATM, ...
               moneyness_modified, c_mkt_calibration, expiries);
end

end
