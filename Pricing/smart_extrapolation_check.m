function R = smart_extrapolation_check(models, params_list, scale_list, ...
        N_sim, M_ref, M_list, dz, forward, K1, K2, df, seed)
% SMART_EXTRAPOLATION_CHECK  Does a smart extrapolation of the simulation CDF
%   outside the FFT range (cf. point 0 and [5]) improve the exotics' MC price
%   on a "not-so-large" grid?  (Project 2, exercise 4.)
%
%   Mechanism. The increment CDF is reconstructed by FFT only on the finite
%   space grid [-zn, zn], zn = dz*N/2. "Smart extrapolation" replaces the
%   truncated tails beyond +-zn with the analytic EXPONENTIAL tails of the
%   model (smart_cdf_extrapolation), so the inverse-CDF sampler can reach
%   beyond the grid instead of piling the tail mass at the edges.
%
%   Caveat made explicit by the sweep. With a single FFT the space range and
%   the Fourier step obey the duality  zn * dx = pi  (dx = 2*pi/(N*dz)). Hence
%   shrinking N at fixed dz shrinks zn AND coarsens dx together: a grid too
%   small to cover the tails is also too coarse to resolve the CF, so the
%   dominant error there is the CF DISCRETISATION, which tail extrapolation
%   cannot fix. We therefore sweep M_sim over a range where the CF is actually
%   resolved and let the numbers show how much (if anything) extrapolation adds.
%
%   For each (model, M_sim) the Call-on-Call is priced on the SAME random
%   stream (rng(seed) before every run) in two configurations:
%     - smart OFF : CDF cleaned but NOT extended (truncated support);
%     - smart ON  : CDF extended with the analytic exponential tails.
%   Both are compared to the large-grid reference (M_ref). The inner-call
%   Lewis-FFT always uses M_ref, so only the increment-sampling grid varies.
%
% INPUTS:
%   models      - cell array of model tags, e.g. {'GL','AB','MA'}
%   params_list - cell array of parameter vectors aligned with `models`
%   scale_list  - cell array of [scale_t1, scale_t2] aligned with `models`
%   N_sim       - number of MC paths
%   M_ref       - reference / inner-call FFT exponent (large, e.g. 16)
%   M_list      - vector of simulation-grid exponents to test (e.g. [11 12 13 14])
%   dz          - z-domain step
%   forward     - F(t0,T2)
%   K1, K2      - compound and inner strikes
%   df          - [B(t0,t1), B(t0,t2)]
%   seed        - RNG seed shared by all runs
%
% OUTPUT:
%   R           - results table (one row per model x M_sim)

    fwd_factor = df(1) / df(2);
    df_t1_t2   = df(2) / df(1);

    Model    = {};
    M_sim    = [];
    N_pts    = [];
    Price_ref      = [];
    Price_smartOFF = [];
    Price_smartON  = [];
    AbsErr_OFF     = [];
    AbsErr_ON      = [];

    for i = 1:numel(models)
        m  = models{i};
        p  = params_list{i};
        sc = scale_list{i};
        [cf_marg, cf_inc] = model_cf_handles(m);

        % --- reference price on the large grid ---------------------------
        rng(seed);
        ft1_ref = sample_ft1(cf_marg, p, sc(1), M_ref, dz, m, N_sim, false);
        p_ref   = price_coc(ft1_ref, cf_inc, p, sc, M_ref, dz, m, ...
                            forward, K1, K2, df, df_t1_t2, fwd_factor);

        for M = M_list(:)'
            % smart OFF: small grid, truncated support
            rng(seed);
            ft1_off = sample_ft1(cf_marg, p, sc(1), M, dz, m, N_sim, false);
            p_off   = price_coc(ft1_off, cf_inc, p, sc, M_ref, dz, m, ...
                                forward, K1, K2, df, df_t1_t2, fwd_factor);
            % smart ON: small grid, analytic tails extended
            rng(seed);
            ft1_on  = sample_ft1(cf_marg, p, sc(1), M, dz, m, N_sim, true);
            p_on    = price_coc(ft1_on, cf_inc, p, sc, M_ref, dz, m, ...
                                forward, K1, K2, df, df_t1_t2, fwd_factor);

            Model{end+1,1}        = m;            %#ok<AGROW>
            M_sim(end+1,1)        = M;            %#ok<AGROW>
            N_pts(end+1,1)        = 2^M;          %#ok<AGROW>
            Price_ref(end+1,1)    = p_ref;        %#ok<AGROW>
            Price_smartOFF(end+1,1) = p_off;      %#ok<AGROW>
            Price_smartON(end+1,1)  = p_on;       %#ok<AGROW>
            AbsErr_OFF(end+1,1)   = abs(p_off - p_ref);  %#ok<AGROW>
            AbsErr_ON(end+1,1)    = abs(p_on  - p_ref);  %#ok<AGROW>
        end
    end

    R = table(Model, M_sim, N_pts, Price_ref, Price_smartOFF, AbsErr_OFF, ...
              Price_smartON, AbsErr_ON);
end

% =========================================================================
%  Helpers
% =========================================================================
function [cf_marg, cf_inc] = model_cf_handles(model)
% Marginal CF at t1 (for sampling ft1) and increment CF t1->t2 (inner call).
    switch upper(model)
        case 'GL', cf_marg = @cf_GL;    cf_inc = @cf_increment_GL;
        case 'AB', cf_marg = @cf_AB;    cf_inc = @cf_increment_AB;
        case 'MA', cf_marg = @cf_MA_IA; cf_inc = @cf_MA_FA;
        otherwise, error('smart_extrapolation_check:badModel', 'Unknown model %s', model);
    end
end

function ft1 = sample_ft1(cf_marg, params, scale_t1, M, dz, model, N_sim, smart)
% Invert the Lewis-FFT digital CDF of the t1 marginal; optionally extend the
% support beyond the FFT range with the analytic exponential tails.
    [cdf_raw, z_raw] = lewis_FFT_digital(cf_marg, M, dz, params, scale_t1, ...
                                         ~smart, model, 1);   % clean within-range when NOT smart
    if smart
        [cdf_use, z_use] = smart_cdf_extrapolation(z_raw, cdf_raw);
    else
        cdf_use = cdf_raw;  z_use = z_raw;
    end
    ft1 = simulate_from_cdf(cdf_use, z_use, 1, N_sim);
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
