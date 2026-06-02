function LA_results = run_ex4(params, market, diagnostics)
%RUN_EX4 Computes prices for Call-on-Call, Put-on-Put, and Chooser options 
%        using Linear Additive Models (GL, MA, AB) and performs sanity checks.

    addpath("Utilities/");
    addpath("Distributions/");
    addpath(genpath('Pricing/'));
    addpath(genpath('Simulation/'));
    addpath(genpath('Calibration/'));

    % Handle default arguments
    if nargin < 3 || isempty(diagnostics)
        diagnostics = false;
    end

    % Standalone mode: if no calibrated inputs are passed, calibrate first.
    if nargin < 2
        [params, market] = calibrate_surface(struct('verbose', false));
    end

    fprintf('=========================================================================\n');
    fprintf('        LINEAR ADDITIVE MODELS (MA / AB / GL) - EXOTIC PRICE ENGINE      \n');
    fprintf('=========================================================================\n\n');

    %% =========================================================================
    % STEP 1-2: UNPACK MARKET DATA AND CALIBRATED PARAMETERS
    % =========================================================================
    fprintf('STEP 1-2: Using pre-calibrated curve, ATM vol and model parameters...\n');
    
    strikes         = market.strikes;
    calls           = market.calls;
    puts            = market.puts;
    expiries        = market.expiries;
    discount_factor = market.discount_factor;
    forward         = market.forward;
    sigma_ATM       = market.sigma_ATM;
    yf              = market.yf;
    nT              = numel(expiries);
    
    fprintf('  -> Market data and ATM volatility loaded for %d maturities.\n\n', nT);
    
    iT1 = 2;  iT2 = 4;
    T1  = yf(iT1);  T2 = yf(iT2);
    fprintf('  Forward-start window: t1 = %.4f y (idx %d), t2 = %.4f y (idx %d)\n\n', T1, iT1, T2, iT2);

    %% =========================================================================
    %  COMMON PARAMETERS (shared across all models)
    %% =========================================================================
    M      = 16;
    dz     = 5e-3;
    N_sim  = 1e7;
    N_grid = 300;
    
    F_t0_t2 = forward(iT2);
    K2      = F_t0_t2;          % ATM strike for vanilla call
    K1      = 1;                % compound option strike
    df      = [discount_factor(iT1), discount_factor(iT2)];
    B_0_t2  = discount_factor(iT2);
    
    % Lemma 2 (Forward.pdf): F(T1,T2) = F(t0,T2) + fwd_factor * f_{T1,T1}. The MC
    % pricers recompute this internally from df; it is needed here only for the
    % GL/AB ATM chooser closed-form decomposition (the t1 leg is on F(T1,T2)).
    fwd_factor = discount_factor(iT1) / discount_factor(iT2);

    % Sanity-check enabling flag: the Chooser closed-form decomposition into
    % two ATM-forward calls is valid ONLY when K2 = F(t0,T2).
    chooser_sanity_enabled = (K2 == F_t0_t2);

    %% =========================================================================
    %  GL MODEL
    %% =========================================================================
    fprintf('--- 1. GENERALIZED LAPLACE (GL) MODEL ---\n');
    params_GL        = params.GL; 
    I0_opt_GL        = I0_GL(params_GL);
    sigma_t_GL       = sigma_ATM / I0_opt_GL;
    scale_factor_GL  = [sigma_t_GL(iT1)*sqrt(yf(iT1)), sigma_t_GL(iT2)*sqrt(yf(iT2))];
    
    rng(1234);
    fprintf('\n  > Computing Call-on-Call (MC)...\n');
    [price_GL, CI_GL, ft1_GL, call_price_t1_GL] = CoC_pricing_MC( ...
        params_GL, scale_factor_GL, N_sim, M, dz, N_grid, F_t0_t2, K1, K2, df, 'GL', diagnostics);
        
    fprintf('\n  > Computing Put-on-Put (MC)...\n');
    [price_PoP_GL, CI_PoP_GL, ft1_PoP_GL, put_price_t1_GL] = PoP_pricing_MC( ...
        params_GL, scale_factor_GL, N_sim, M, dz, N_grid, F_t0_t2, K1, K2, df, 'GL', diagnostics);
        
    fprintf('\n  > Computing Chooser (MC)...\n');
    [price_Ch_GL, CI_Ch_GL, ft1_Ch_GL, call_price_Ch_t1_GL] = Chooser_pricing_MC( ...
        params_GL, scale_factor_GL, N_sim, M, dz, N_grid, F_t0_t2, K2, df, 'GL', diagnostics);
        
    fprintf('  > Computing Sanity Checks (K1=0 CoC & Chooser Analytic)...\n\n');
    rng(1234);
    [price_CoC_K1zero_GL, ~] = CoC_pricing_MC( ...
        params_GL, scale_factor_GL, N_sim, M, dz, N_grid, F_t0_t2, 0, K2, df, 'GL', false);
    [price_ATM_GL] = call_ATM_vanilla(params_GL, scale_factor_GL(2), B_0_t2, 'GL');
    
    if chooser_sanity_enabled
        % t1 leg is the ATM option on F(T1,T2) = F(t0,T2) + fwd_factor*f_{T1,T1},
        % hence the scale carries fwd_factor (Lemma 2). call_ATM_vanilla is linear
        % in the scale, so this equals fwd_factor * ATM-call(T1).
        [call_ATM_T2_GL] = call_ATM_vanilla(params_GL, scale_factor_GL(2), B_0_t2, 'GL');
        [call_ATM_T1_GL] = call_ATM_vanilla(params_GL, fwd_factor*scale_factor_GL(1), B_0_t2, 'GL');
        price_Ch_GL_analytic    = call_ATM_T2_GL + call_ATM_T1_GL;
    else
        price_Ch_GL_analytic = NaN;
    end

    %% =========================================================================
    %  MA MODEL
    %% =========================================================================
    fprintf('--- 2. MINIMAL ADDITIVE (MA) MODEL ---\n');
    params_MA        = params.MA; 
    I0_opt_MA        = I0_MA(params_MA);
    sigma_t_MA       = sigma_ATM / I0_opt_MA;
    scale_factor_MA  = [sigma_t_MA(iT1)*sqrt(yf(iT1)), sigma_t_MA(iT2)*sqrt(yf(iT2))];
    
    rng(1234);
    fprintf('\n  > Computing Call-on-Call (MC, Semi-Analytic, Fully Analytic)...\n');
    [price_CoC_MA_MC, CI_MA, ft1_MA, call_price_t1_MA] = CoC_pricing_MC( ...
        params_MA, scale_factor_MA, N_sim, M, dz, N_grid, F_t0_t2, K1, K2, df, 'MA', diagnostics);
    price_CoC_MA_analytic = CoC_pricing_analytical(params_MA, scale_factor_MA, F_t0_t2, K1, K2, df);
    price_CoC_MA_FULLanalytic = CoC_pricing_FULLanalytical(params_MA, scale_factor_MA, F_t0_t2, K1, K2, df);
    diff_CoC_MA_bps = (price_CoC_MA_FULLanalytic - price_CoC_MA_analytic) * 10000;
    
    if diagnostics
        fprintf('\n=================================================================\n');
        fprintf('          CoC MA: SEMI-ANALYTIC VS FULLY ANALYTIC                \n');
        fprintf('=================================================================\n');
        fprintf('  Semi-Analytic Price (quadgk) : %12.8f\n', price_CoC_MA_analytic);
        fprintf('  Fully Analytic Price (Exact) : %12.8f\n', price_CoC_MA_FULLanalytic);
        fprintf('  Difference                   : %12.4f bps\n', diff_CoC_MA_bps);
        fprintf('=================================================================\n\n');
    end
    
    fprintf('\n  > Computing Put-on-Put (MC & Analytical)...\n');
    [price_PoP_MA_MC, CI_PoP_MA, ft1_PoP_MA, put_price_t1_MA] = PoP_pricing_MC( ...
        params_MA, scale_factor_MA, N_sim, M, dz, N_grid, F_t0_t2, K1, K2, df, 'MA', diagnostics);
    price_PoP_MA_analytic = PoP_pricing_analytical(params_MA, scale_factor_MA, F_t0_t2, K1, K2, df, diagnostics);
    
    fprintf('\n  > Computing Chooser (MC & Analytical)...\n');
    [price_Ch_MA_MC, CI_Ch_MA, ft1_Ch_MA, call_price_Ch_t1_MA] = Chooser_pricing_MC( ...
        params_MA, scale_factor_MA, N_sim, M, dz, N_grid, F_t0_t2, K2, df, 'MA', diagnostics);
    price_Ch_MA_analytic = Chooser_pricing_analytic(params_MA, scale_factor_MA, df, F_t0_t2, K2);
    
    if diagnostics
        K_eval = linspace(-10, 10, 20);
        fprintf('\n>>> 1. FINITE ACTIVITY INCREMENT (Delta f)\n');
        fprintf('Comparing analytical formula vs FFT Lewis inversion...\n');
        comparison_function_increments_MA(params_MA, scale_factor_MA, df, M, dz, K_eval, F_t0_t2);
        
        fprintf('\n>>> 2. INFINITE ACTIVITY MARGINAL (f_t)\n');
        fprintf('Comparing analytical formula vs FFT Lewis inversion...\n');
        comparison_function_marginals_MA(params_MA, scale_factor_MA, df, M, dz, K_eval);
    end
    fprintf('\n');

    %% =========================================================================
    %  AB MODEL
    %% =========================================================================
    fprintf('--- 3. ADDITIVE BACHELIER (AB) MODEL ---\n');
    params_AB        = params.AB; 
    I0_opt_AB        = I0_AB(0, params_AB);
    sigma_t_AB       = sigma_ATM / I0_opt_AB;
    scale_factor_AB  = [sigma_t_AB(iT1)*sqrt(yf(iT1)), sigma_t_AB(iT2)*sqrt(yf(iT2))];
    
    rng(1234);
    fprintf('\n  > Computing Call-on-Call (MC)...\n');
    [price_AB, CI_AB, ft1_AB, call_price_t1_AB] = CoC_pricing_MC( ...
        params_AB, scale_factor_AB, N_sim, M, dz, N_grid, F_t0_t2, K1, K2, df, 'AB', diagnostics);
        
    fprintf('\n  > Computing Put-on-Put (MC)...\n');
    [price_PoP_AB, CI_PoP_AB, ft1_PoP_AB, put_price_t1_AB] = PoP_pricing_MC( ...
        params_AB, scale_factor_AB, N_sim, M, dz, N_grid, F_t0_t2, K1, K2, df, 'AB', diagnostics);
        
    fprintf('\n  > Computing Chooser (MC)...\n');
    [price_Ch_AB, CI_Ch_AB, ft1_Ch_AB, call_price_Ch_t1_AB] = Chooser_pricing_MC( ...
        params_AB, scale_factor_AB, N_sim, M, dz, N_grid, F_t0_t2, K2, df, 'AB', diagnostics);
        
    fprintf('\n  > Computing Sanity Checks (K1=0 CoC & Chooser Analytic)...\n\n');
    rng(1234);
    [price_CoC_K1zero_AB, ~] = CoC_pricing_MC( ...
        params_AB, scale_factor_AB, N_sim, M, dz, N_grid, F_t0_t2, 0, K2, df, 'AB', false);
    [price_ATM_AB] = call_ATM_vanilla(params_AB, scale_factor_AB(2), B_0_t2, 'AB');
    
    if chooser_sanity_enabled
        % t1 leg on F(T1,T2): scale carries fwd_factor (Lemma 2), see GL block above.
        [call_ATM_T2_AB] = call_ATM_vanilla(params_AB, scale_factor_AB(2), B_0_t2, 'AB');
        [call_ATM_T1_AB] = call_ATM_vanilla(params_AB, fwd_factor*scale_factor_AB(1), B_0_t2, 'AB');
        price_Ch_AB_analytic    = call_ATM_T2_AB + call_ATM_T1_AB;
    else
        price_Ch_AB_analytic = NaN;
    end

    %% =========================================================================
    %  RESULTS TABLE — PRICES
    %% =========================================================================
    models = {'GL'; 'MA'; 'AB'};
    
    % --- Call-on-Call Data ---
    prices_CoC_MC       = [price_GL; price_CoC_MA_MC; price_AB];
    prices_CoC_Analytic = [NaN;      price_CoC_MA_FULLanalytic; NaN];
    IC_CoC              = [CI_GL, CI_MA, CI_AB]'; 
    CI_CoC_str          = compose('[%.4f, %.4f]', IC_CoC(:,1), IC_CoC(:,2));
    Err_CoC_bps         = (prices_CoC_Analytic - prices_CoC_MC) * 10000;
    
    % --- Put-on-Put Data ---
    prices_PoP_MC       = [price_PoP_GL; price_PoP_MA_MC; price_PoP_AB];
    prices_PoP_Analytic = [NaN;          price_PoP_MA_analytic; NaN]; 
    IC_PoP              = [CI_PoP_GL, CI_PoP_MA, CI_PoP_AB]'; 
    CI_PoP_str          = compose('[%.4f, %.4f]', IC_PoP(:,1), IC_PoP(:,2));
    Err_PoP_bps         = (prices_PoP_Analytic - prices_PoP_MC) * 10000;  
    
    % --- Chooser Data ---
    prices_Ch_MC       = [price_Ch_GL;          price_Ch_MA_MC;       price_Ch_AB];
    prices_Ch_Analytic = [price_Ch_GL_analytic; price_Ch_MA_analytic; price_Ch_AB_analytic];
    IC_Ch              = [CI_Ch_GL, CI_Ch_MA, CI_Ch_AB]'; 
    CI_Ch_str          = compose('[%.4f, %.4f]', IC_Ch(:,1), IC_Ch(:,2));
    Err_Ch_bps         = (prices_Ch_Analytic - prices_Ch_MC) * 10000;
    
    % --- Create and Display Table ---
    T = table(models, ...
              prices_CoC_MC, prices_CoC_Analytic, Err_CoC_bps, CI_CoC_str, ...
              prices_PoP_MC, prices_PoP_Analytic, Err_PoP_bps, CI_PoP_str, ...
              prices_Ch_MC, prices_Ch_Analytic, Err_Ch_bps, CI_Ch_str, ...
        'VariableNames', {'Model', ...
                          'CoC_MC', 'CoC_Analytic', 'CoC_Err_bps', 'CoC_95_CI', ...
                          'PoP_MC', 'PoP_Analytic', 'PoP_Err_bps', 'PoP_95_CI', ...
                          'Chooser_MC', 'Chooser_Analytic', 'Chooser_Err_bps', 'Chooser_95_CI'});
                          
    disp(' ');
    disp('=================================================================================================================================');
    disp('                                                        PRICING RESULTS                                                          ');
    disp('=================================================================================================================================');
    disp(T);

    %% =========================================================================
    %  SANITY CHECKS SUMMARY 
    %
    %  Test 1: CoC with K1 = 0 must collapse to a vanilla ATM-forward call at T2.
    %          Reference: call_ATM_vanilla at T2.
    %          Valid for any K2 (the K1=0 collapse does not require ATM strike).
    %
    %  Test 2: Chooser = ATM-call(T2) + ATM-call(T1).
    %          Reference: sum of two call_ATM_vanilla evaluations.
    %          Valid ONLY when K2 = F(t0,T2) (ATM). Outside this case the
    %          closed-form does not hold and the check is skipped.
    %
    %  MA is excluded from these checks: it has its own closed-form analytic
    %  pricer already validated against MC in the table above.
    %% =========================================================================
    err_CoC_GL_bps = (price_ATM_GL - price_CoC_K1zero_GL) * 10000;
    err_CoC_AB_bps = (price_ATM_AB - price_CoC_K1zero_AB) * 10000;
    
    if chooser_sanity_enabled
        err_Ch_GL_bps = (price_Ch_GL_analytic - price_Ch_GL) * 10000;
        err_Ch_AB_bps = (price_Ch_AB_analytic - price_Ch_AB) * 10000;
    else
        err_Ch_GL_bps = NaN;
        err_Ch_AB_bps = NaN;
    end
    
    sanity_models = {'GL'; 'AB'};
    sanity_CoC_MC = [price_CoC_K1zero_GL; price_CoC_K1zero_AB];
    sanity_CoC_Analytic = [price_ATM_GL;  price_ATM_AB];
    sanity_CoC_Er_bps = [err_CoC_GL_bps;  err_CoC_AB_bps];
    
    sanity_Ch_MC  = [price_Ch_GL;           price_Ch_AB];
    sanity_Ch_Analytic  = [price_Ch_GL_analytic;  price_Ch_AB_analytic];
    sanity_Ch_Er_bps  = [err_Ch_GL_bps;     err_Ch_AB_bps];
    
    S = table(sanity_models, sanity_CoC_MC, sanity_CoC_Analytic, sanity_CoC_Er_bps, ...
                             sanity_Ch_MC,  sanity_Ch_Analytic,  sanity_Ch_Er_bps, ...
        'VariableNames', {'Model', ...
                          'CoC_K1zero_MC', 'CoC_VanillaATM_Analytic', 'CoC_Err_bps', ...
                          'Chooser_MC',    'Chooser_SR_Analytic',     'Chooser_Err_bps'});
                          
    disp(' ');
    disp('=============================================================================================================');
    disp('                                       SANITY CHECKS                                                         ');
    disp('  CoC: K1=0 collapse vs vanilla ATM call at T2   (always valid)                                              ');
    if chooser_sanity_enabled
        disp('  Chooser: ATM-call(T2) + ATM-call(T1) vs MC      (valid because K2 = F(t0,T2))                              ');
    else
        disp('  Chooser: skipped (closed-form requires K2 = F(t0,T2), current K2 differs).                                 ');
    end
    disp('=============================================================================================================');
    disp(S);
    
    %% =========================================================================
    %  SMART EXTRAPOLATION OF THE SIMULATION CDF — TRUNCATION TEST
    %
    %  Ex. 4 hint (cf. point 0 and [5]): a "smart extrapolation" replaces the CDF
    %  tails truncated by the finite FFT range with the model's analytic exponential
    %  tails. To isolate the truncation effect from the (unrelated) CF discretisation
    %  error, we start from the ACCURATE increment CDF (large grid M_ref) and
    %  ARTIFICIALLY truncate it to a symmetric window [-L, L], L = c * std(ft1). For
    %  each truncation level the CoC is priced with smart OFF (tails lost) vs ON
    %  (analytic tails re-attached) on the same random stream, vs the untruncated
    %  reference. AbsErr_OFF should grow as L shrinks, AbsErr_ON should stay small:
    %  this proves the extrapolation works. In the real FFT grid such truncation
    %  never arises (a CF-resolving grid spans >> 10 std), so for these light-tailed
    %  models the extrapolation is in practice unnecessary.
    %% =========================================================================
    c_list    = [1, 1.5, 2, 3];   % truncation half-widths in std units (L = c*std)
    N_sim_ext = 1e6;              % shared-seed runs -> MC noise cancels
    models_ext = {'GL'; 'AB'; 'MA'};
    params_ext = {params_GL; params_AB; params_MA};
    scale_ext  = {scale_factor_GL; scale_factor_AB; scale_factor_MA};
    
    R_ext = smart_extrapolation_check(models_ext, params_ext, scale_ext, ...
                N_sim_ext, M, c_list, dz, F_t0_t2, K1, K2, df, 1234);
                
    disp(' ');
    disp('=============================================================================================');
    disp('               SMART EXTRAPOLATION OF THE SIMULATION CDF — TRUNCATION TEST                    ');
    disp('  Accurate CDF truncated to [-L, L] (L = c*std): smart OFF (tails lost) vs ON (analytic       ');
    disp('  tails) vs untruncated reference. TailCut = discarded mass; AbsErr = |price - reference|.    ');
    disp('=============================================================================================');
    fprintf('  Reference grid M = %d (N = %d)  |  N_sim = %.0e\n\n', M, 2^M, N_sim_ext);
    disp(R_ext);
    
    % Pack tables into the output structure
    LA_results.Pricing            = T;
    LA_results.SanityChecks       = S;
    LA_results.SmartExtrapolation = R_ext;
end
