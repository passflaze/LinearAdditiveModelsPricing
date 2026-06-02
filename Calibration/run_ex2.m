function [params, market] = run_ex2(opts)
% CALIBRATE_SURFACE Calibrates the volatility surface to the 3 Linear Additive
% Models (Additive Bachelier, Minimal Additive, Generalized Laplace) and
% returns the calibrated parameters.
%
% This is the "function" version of the ex2.m script: it performs curve
% bootstrapping, ATM vol calibration, surface filtering, and calibration
% of the 3 models, returning the parameters in a reusable struct (e.g., for
% run_ex3.m) without having to repeat the code.
%
% INPUT (opts, optional struct -- all fields have defaults):
%   opts.callpath    : calls file path            (default "Data/datacalls")
%   opts.putpath     : puts file path             (default "Data/dataputs")
%   opts.expiryFile  : expiries file path         (default "Data/Expiries_Futures.txt")
%   opts.valueDate   : value date (datetime)      (default 2020-06-02)
%   opts.x_min       : dollar-moneyness band min  (default -30)
%   opts.x_max       : dollar-moneyness band max  (default  30)
%   opts.M           : FFT exponent for GL        (default 15)
%   opts.dz          : FFT grid step for GL       (default 2.5e-3)
%   opts.verbose     : true to print the report   (default true)
%   opts.plot        : true to plot distributions (default false)
%
% OUTPUT:
%   params : struct with calibrated parameters; each field is a COLUMN VECTOR
%            according to the params(1)/params(2) convention used throughout
%            the project (cf_*, pdf_*, price_*, lewis_FFT_*).
%       .AB   -> [k;     eta ]
%       .MA   -> [alpha; beta]
%       .GL   -> [alpha; beta]
%   market : struct with market and support data (column vectors)
%       .strikes, .calls, .puts, .expiries, .valueDate
%       .discount_factor, .forward, .R2
%       .c_ATM, .sigma_ATM, .yf
%       .moneyness_modified, .c_mkt_calibration

addpath("Calibration/Calibration_MA/");
addpath("Calibration/Calibration_AB/");
addpath("Calibration/Calibration_GL/");

% =========================================================================
% INPUT HANDLING / DEFAULTS
% =========================================================================
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
if ~isfield(opts, 'verbose'),    opts.verbose    = true;                        end
if ~isfield(opts, 'plot'),       opts.plot       = false;                       end

vb = @(varargin) verbose_print(opts.verbose, varargin{:});
vb('=========================================================================\n');
vb('             VOLATILITY SURFACE CALIBRATION ENGINE                       \n');
vb('=========================================================================\n\n');

%% =========================================================================
% STEP 1: DATA LOADING AND BOOTSTRAPPING
% =========================================================================
vb('STEP 1: Loading data and bootstrapping curve...\n');

% Read market options data
[strikes, calls, puts, expiries] = readData(opts.callpath, opts.putpath, ...
                                             opts.valueDate, opts.expiryFile);

% Bootstrap synthetic discount factors and forwards
nT = numel(expiries);
discount_factor = zeros(nT, 1);
forward         = zeros(nT, 1);
R2              = zeros(nT, 1);

for k = 1:nT
    [discount_factor(k), forward(k), R2(k)] = bootstrap(puts(k,:), calls(k,:), strikes);
end

% Report Bootstrap Results
vb("\n--- Bootstrap Results (Value Date: %s) ---\n", string(opts.valueDate, "yyyy-MM-dd"));
vb("%-12s  %10s  %10s  %8s\n", "Expiry", "D(t,T)", "F(t,T)", "R^2");
for k = 1:nT
    vb("%-12s  %10.6f  %10.4f  %8.4f\n", ...
        string(expiries(k), "yyyy-MM-dd"), ...
        discount_factor(k), forward(k), R2(k));
end
vb('-------------------------------------------------------------------------\n\n');

% Calculate the standard year fraction and zero rates
yf = yearfrac(opts.valueDate, expiries, 3);
zero_rates = -log(discount_factor) ./ yf;

%% =========================================================================
% STEP 2: ATM VOLATILITY CALIBRATION & MONEYNESS GENERATION
% =========================================================================
vb('STEP 2: Calibrating ATM Volatility and filtering surface...\n');
c_ATM = zeros(length(forward), 1);
for i = 1:length(forward)
    current_calls = calls(i, :);
    current_puts  = puts(i, :);
    c_ATM(i) = callATM(current_calls, current_puts, strikes, forward(i), discount_factor(i));
end

sigma_ATM = sigmaATM(c_ATM, discount_factor, yf, expiries);
check_term_structure(sigma_ATM, yf, expiries);

x_min = opts.x_min;
x_max = opts.x_max;
[moneyness_modified, c_mkt_calibration] = moneyness_generator( ...
    forward, strikes, calls, puts, sigma_ATM, yf, discount_factor, ...
    x_min, x_max);
vb('  -> Surface filtered on dollar moneyness x in [%g, %g] $ and prepared for optimization.\n\n', ...
    x_min, x_max);

%% =========================================================================
% STEP 3: ADDITIVE BACHELIER (AB) MODEL CALIBRATION
% =========================================================================
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

%% =========================================================================
% STEP 4: MINIMAL ADDITIVE (MA) MODEL CALIBRATION
% =========================================================================
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

check_skew_MA(alpha_MA, beta_MA, discount_factor, yf, ...
    sigma_ATM, moneyness_modified, c_mkt_calibration, expiries);

%% =========================================================================
% STEP 5: GENERALIZED LAPLACE (GL) MODEL CALIBRATION
% =========================================================================
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

%% =========================================================================
% STEP 6: DIAGNOSTICS & PACK OUTPUTS
% =========================================================================
if opts.verbose
    print_diagnostics(k_AB, eta_AB, fval_AB, alpha_MA, beta_MA, fval_MA, ...
        alpha_GL, beta_GL, fval_GL, M, dz, discount_factor, yf, sigma_ATM, ...
        moneyness_modified, c_mkt_calibration, expiries, nT, opts.plot);
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
    market.zero_rates         = zero_rates; % Added for completeness
end

%% =========================================================================
% STEP 7: PLOTTING SUITE (Executed only if opts.plot == true)
% =========================================================================
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

    % --- 7.2: Volatility Smiles Comparison ---
    sigma_MA = sigma_ATM ./ I0_MA(params.MA);
    sigma_AB = sigma_ATM ./ I0_AB(0, params.AB);
    sigma_GL = sigma_ATM ./ I0_GL(params.GL);
    
    figure('Name', 'Volatility Smiles Comparison', 'Color', 'w');
    hold on;
    plot(expiries, sigma_MA, 'r-', 'LineWidth', 1.5, 'DisplayName', 'MA Volatility');
    plot(expiries, sigma_AB, 'b-', 'LineWidth', 1.5, 'DisplayName', 'AB Volatility');
    plot(expiries, sigma_GL, 'g-', 'LineWidth', 1.5, 'DisplayName', 'GL Volatility');
    plot(expiries, sigma_ATM, 'k--', 'LineWidth', 1.5, 'DisplayName', 'ATM Volatility');
    xlabel('Maturity');
    ylabel('Implied Volatility');
    title('Volatility Smiles Comparison vs ATM Volatility');
    legend('Location', 'best');
    grid on;
    hold off;

    % --- 7.3: Distributions and Implied Volatility Smile (Bachelier) ---
    plot_distributions(eta_AB, k_AB, alpha_MA, beta_MA, alpha_GL, beta_GL);
    
    plot_iv_AB(k_AB, eta_AB, discount_factor, yf, sigma_ATM, ...
               moneyness_modified, c_mkt_calibration, expiries);
end

end
