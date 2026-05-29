function results = run_ex3_simulation(params, opts)
%RUN_EX3_SIMULATION  Exercise 3 driver: forward-start option pricing under the
%   three Linear Additive models (MA, AB, GL), comparing analytic vs Monte Carlo.
%
%   results = RUN_EX3_SIMULATION(params)
%   results = RUN_EX3_SIMULATION(params, opts)
%
%   It takes the already-CALIBRATED shape parameters as input (decoupling ex3
%   from the ex2 calibration); the market curve and ATM vols are rebuilt from
%   data, and each model scale sigma_t = sigma_ATM/I0 is derived internally.
%
%     - MA: analytic vs MC across a set of proportional strikes K2 (exact
%           closed form available for every strike).
%     - AB, GL: analytic vs MC at the ATM forward-start (K2 = 1), where the
%           payoff reduces to max(W,0) and the analytic price is the survival
%           integral E[max(W,0)] = int_0^inf (1 - F_W) dx. The two twin models
%           run through one identical, unified engine.
%
%   INPUT  params : struct of calibrated shape parameters
%            params.MA.alpha, params.MA.beta
%            params.GL.alpha, params.GL.beta
%            params.AB.k    , params.AB.eta
%   INPUT  opts   : optional struct (defaults applied)
%            .iT1 (2), .iT2 (4)   reset / maturity indices
%            .N_sim (1e6), .seed (2), .n_strikes (20)  (MA strike set)
%            .doPlots (true)      MA analytic-vs-MC curve
%            .callpath/.putpath/.expiryFile/.valueDate  (Data/ defaults)

    % ---- paths & options ------------------------------------------------
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(here, fullfile(root,'Utilities'));

    if nargin < 2, opts = struct(); end
    opts = fill_defaults(opts, struct( ...
        'callpath',   fullfile(root,'Data','datacalls'), ...
        'putpath',    fullfile(root,'Data','dataputs'), ...
        'expiryFile', fullfile(root,'Data','Expiries_Futures.txt'), ...
        'valueDate',  datetime(2020,6,2), ...
        'iT1', 2, 'iT2', 4, 'N_sim', 1e6, 'seed', 2, ...
        'n_strikes', 20, 'doPlots', true));

    % ---- market curve (deterministic) -----------------------------------
    [strikes, calls, puts, expiries] = readData(opts.callpath, opts.putpath, ...
                                                opts.valueDate, opts.expiryFile);
    nT = numel(expiries);
    [discount_factor, forward, call_atm] = deal(zeros(nT,1));
    for k = 1:nT
        [discount_factor(k), forward(k), ~] = bootstrap(puts(k,:), calls(k,:), strikes);
        call_atm(k) = callATM(calls(k,:), puts(k,:), strikes, forward(k), discount_factor(k));
    end
    yf        = yearfrac(opts.valueDate, expiries, 3);
    sigma_ATM = sigmaATM(call_atm, discount_factor, yf, expiries);

    iT1 = opts.iT1;  iT2 = opts.iT2;
    T1  = yf(iT1);   T2 = yf(iT2);
    B1  = discount_factor(iT1);  B2 = discount_factor(iT2);
    F2  = forward(iT2);

    fprintf('Forward-start window: T1 = %.3f y, T2 = %.3f y  (F = %.4f, B = %.6f)\n', ...
            T1, T2, F2, B2);

    % =====================================================================
    % MINIMAL ADDITIVE (MA): analytic vs MC across strikes
    % =====================================================================
    rng(opts.seed);
    K2 = sort([1, 0.8 + 0.4*rand(1, opts.n_strikes)]);     % proportional strikes

    aMA = params.MA.alpha;  bMA = params.MA.beta;
    gammaMA = 1/aMA - 1/bMA;
    C_MA    = 1 / (1/aMA + 1/bMA);
    I0MA    = I0_MA(gammaMA, C_MA, aMA, bMA);
    sigmat  = (sigma_ATM([iT1 iT2]) / I0MA) .* yf([iT1 iT2]);

    MA_an = pricing_fwd_start_analytic(aMA, bMA, sigmat, B2, K2, F2);
    rng(opts.seed);
    [MA_mc, MA_ci] = pricing_fwd_start_MC(F2, K2, B2, opts.N_sim, 16, 5e-3, ...
                                          sigmat, aMA, bMA);

    fprintf('\n=== MA: forward-start analytic vs Monte Carlo ===\n');
    disp(table(K2(:), MA_an(:), MA_mc(:), MA_ci(1,:)', MA_ci(2,:)', (MA_mc-MA_an)'*1e4, ...
        'VariableNames', {'K2','Analytic','MC','CI_low','CI_high','Diff_bps'}));

    if opts.doPlots
        [Ks, o] = sort(K2);
        figure('Name','MA forward-start'); hold on; grid on;
        plot(Ks, MA_an(o), 'b-', 'LineWidth', 2, 'DisplayName','Analytic');
        plot(Ks, MA_mc(o), 'ro', 'MarkerSize', 6, 'DisplayName','Monte Carlo');
        xlabel('Proportional strike K_2'); ylabel('Discounted price');
        title('MA forward-start: analytic vs MC'); legend('Location','best'); hold off;
    end

    % =====================================================================
    % ADDITIVE BACHELIER (AB) & GENERALIZED LOGISTIC (GL): K2 = 1 only
    %   same unified engine, one identical loop
    % =====================================================================
    twin_models = {'AB','GL'};
    nM   = numel(twin_models);
    twin = struct();
    rows = strings(nM,1);
    [an, mc, lo, hi] = deal(zeros(nM,1));
    for m = 1:nM
        name = twin_models{m};
        spec = la_model_spec(name, params.(name));
        s1   = sigma_ATM(iT1) / spec.I0;
        s2   = sigma_ATM(iT2) / spec.I0;

        rng(opts.seed);
        [p_mc, ci, p_an] = price_fwd_start_LA(spec, T1, T2, s1, s2, F2, B1, B2, 1, opts.N_sim);

        rows(m) = name;
        an(m) = p_an;  mc(m) = p_mc;  lo(m) = ci(1);  hi(m) = ci(2);
        twin.(name) = struct('analytic', p_an, 'MC', p_mc, 'CI', ci(:)');
    end

    fprintf('\n=== AB & GL: forward-start analytic vs Monte Carlo (K2 = 1) ===\n');
    disp(table(rows, an, mc, lo, hi, (mc-an)*1e4, ...
        'VariableNames', {'Model','Analytic','MC','CI_low','CI_high','Diff_bps'}));

    % ---- results --------------------------------------------------------
    results = struct();
    results.window = struct('iT1',iT1,'iT2',iT2,'T1',T1,'T2',T2,'forward',F2,'B',B2);
    results.K2 = K2(:);
    results.MA = struct('analytic', MA_an(:), 'MC', MA_mc(:), 'CI', MA_ci');
    results.AB = twin.AB;
    results.GL = twin.GL;
end

% -------------------------------------------------------------------------
function s = fill_defaults(s, def)
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(s, f{i}) || isempty(s.(f{i})), s.(f{i}) = def.(f{i}); end
    end
end
