function R = smart_extrapolation_check(models, params_list, scale_list, ...
        N_sim, M_ref, c_list, dz, forward, K1, K2, df, seed)
% SMART_EXTRAPOLATION_CHECK  Truncation test for the smart extrapolation of the
%   simulation CDF outside the FFT range (cf. point 0 and [5]) — Project 2, ex 4.
%
%   Goal. Isolate the effect of a TRUNCATED CDF support from the (unrelated) CF
%   discretisation error. Shrinking the FFT grid at fixed dz would entangle the
%   two (duality zn*dx = pi). Instead we start from the ACCURATE increment CDF
%   computed on a large grid (M_ref, CF fully resolved) and ARTIFICIALLY truncate
%   it to a symmetric window [-L, L], keeping the correct shape inside. L is set
%   as a multiple c of the increment's standard deviation, L = c * std(ft1).
%
%   For each model and each truncation level the Call-on-Call is priced on the
%   SAME random stream (rng(seed) before every run, variance-free comparison):
%     - reference : sample from the full (untruncated) accurate CDF;
%     - smart OFF : sample from the truncated CDF (mass piled at +-L, tails lost);
%     - smart ON  : truncated CDF re-extended beyond +-L with the analytic
%                   exponential tails (smart_cdf_extrapolation).
%   The inner-call Lewis-FFT always uses M_ref. AbsErr = |price - reference|.
%
%   Expected reading: AbsErr_OFF grows as L shrinks (tail mass lost), while
%   AbsErr_ON stays small (the analytic tails recover it) — provided L is far
%   enough out for the local decay rate to be in the exponential regime. This
%   demonstrates that the extrapolation works; in the real FFT grid such a
%   truncation never arises (a CF-resolving grid already spans >> 10 std), so
%   for these light-tailed models the extrapolation is in practice unnecessary.
%
% INPUTS:
%   models      - cell array of model tags, e.g. {'GL','AB','MA'}
%   params_list - cell array of parameter vectors aligned with `models`
%   scale_list  - cell array of [scale_t1, scale_t2] aligned with `models`
%   N_sim       - number of MC paths
%   M_ref       - reference / inner-call FFT exponent (large, e.g. 16)
%   c_list      - vector of truncation half-widths in std units (e.g. [1 1.5 2 3])
%   dz          - z-domain step
%   forward     - F(t0,T2)
%   K1, K2      - compound and inner strikes
%   df          - [B(t0,t1), B(t0,t2)]
%   seed        - RNG seed shared by all runs
%
% OUTPUT:
%   R           - results table (one row per model x truncation level)

    fwd_factor = df(1) / df(2);
    df_t1_t2   = df(2) / df(1);

    Model       = {};
    L_over_std  = [];
    TailCut_pct = [];
    Price_ref      = [];
    Price_smartOFF = [];
    AbsErr_OFF     = [];
    Price_smartON  = [];
    AbsErr_ON      = [];

    for i = 1:numel(models)
        m  = models{i};
        p  = params_list{i};
        sc = scale_list{i};
        [cf_marg, cf_inc] = model_cf_handles(m);

        % --- accurate CDF on the large grid (CF fully resolved) ----------
        [cdf_acc, z_acc] = lewis_FFT_digital(cf_marg, M_ref, dz, p, sc(1), ...
                                             true, m, 1);
        cdf_acc = cdf_acc(:);  z_acc = z_acc(:);

        % --- reference price + increment std -----------------------------
        rng(seed);
        ft1_ref = simulate_from_cdf(cdf_acc, z_acc, 1, N_sim);
        p_ref   = price_coc(ft1_ref, cf_inc, p, sc, M_ref, dz, m, ...
                            forward, K1, K2, df, df_t1_t2, fwd_factor);
        std_ft1 = std(ft1_ref);

        for c = c_list(:)'
            L    = c * std_ft1;
            mask = abs(z_acc) <= L;
            z_t   = z_acc(mask);
            cdf_t = cdf_acc(mask);

            % truncated tail mass actually discarded (interpolated, robust)
            tail = interp1(z_acc, cdf_acc, -L, 'pchip') + ...
                   (1 - interp1(z_acc, cdf_acc, L, 'pchip'));

            % smart OFF: hard truncation (mass piled at the window edges)
            cdf_off = cdf_t;  cdf_off(1) = 0;  cdf_off(end) = 1;
            rng(seed);
            ft1_off = simulate_from_cdf(cdf_off, z_t, 1, N_sim);
            p_off   = price_coc(ft1_off, cf_inc, p, sc, M_ref, dz, m, ...
                                forward, K1, K2, df, df_t1_t2, fwd_factor);

            % smart ON: re-extend the analytic exponential tails beyond +-L
            [cdf_on, z_on] = smart_cdf_extrapolation(z_t, cdf_t);
            rng(seed);
            ft1_on = simulate_from_cdf(cdf_on, z_on, 1, N_sim);
            p_on   = price_coc(ft1_on, cf_inc, p, sc, M_ref, dz, m, ...
                               forward, K1, K2, df, df_t1_t2, fwd_factor);

            Model{end+1,1}          = m;                 %#ok<AGROW>
            L_over_std(end+1,1)     = c;                 %#ok<AGROW>
            TailCut_pct(end+1,1)    = 100 * tail;        %#ok<AGROW>
            Price_ref(end+1,1)      = p_ref;             %#ok<AGROW>
            Price_smartOFF(end+1,1) = p_off;             %#ok<AGROW>
            AbsErr_OFF(end+1,1)     = abs(p_off - p_ref);%#ok<AGROW>
            Price_smartON(end+1,1)  = p_on;              %#ok<AGROW>
            AbsErr_ON(end+1,1)      = abs(p_on  - p_ref);%#ok<AGROW>
        end
    end

    R = table(Model, L_over_std, TailCut_pct, Price_ref, ...
              Price_smartOFF, AbsErr_OFF, Price_smartON, AbsErr_ON);
end

% --- Local helpers ---
function [cf_marg, cf_inc] = model_cf_handles(model)
% Marginal CF at t1 (for sampling ft1) and increment CF t1->t2 (inner call).
    switch upper(model)
        case 'GL', cf_marg = @cf_GL;    cf_inc = @cf_increment_GL;
        case 'AB', cf_marg = @cf_AB;    cf_inc = @cf_increment_AB;
        case 'MA', cf_marg = @cf_MA_IA; cf_inc = @cf_MA_FA;
        otherwise, error('smart_extrapolation_check:badModel', 'Unknown model %s', model);
    end
end

function price = price_coc(ft1, cf_inc, params, scale_factor, M, dz, model, ...
                           forward, K1, K2, df, df_t1_t2, fwd_factor)
% Semi-analytic CoC value for a given vector of simulated increments ft1.
    F_t1_T2       = forward + fwd_factor * ft1;        % Lemma 2
    strikes       = K2 - F_t1_T2;
    call_price_t1 = df_t1_t2 * lewis_FFT_call(cf_inc, M, dz, params, ...
                        scale_factor, strikes, 1, model, fwd_factor);
    payoff        = max(call_price_t1 - K1, 0);
    price         = mean(df(1) * payoff);
end
