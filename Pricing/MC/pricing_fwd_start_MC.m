function [price, CI, diag, sigma] = pricing_fwd_start_MC(model, params, sigma_T1, sigma_T2, ...
        N_sim, M, dz, forward, B_0_t1, B_0_t2, strike, opts)
% PRICING_FWD_START_MC  Monte Carlo price of the forward-start option
%   payoff = [ S(T2) - K2 * F(T1,T2) ]_+
% under a Linear Additive model. UNIFIED entry point for the three models
% via a switch on `model` ('MA' | 'AB' | 'GL').
%
% Strategy (shared by all models, exploits additivity and the Lemma 2 forward
% rescaling, Forward.pdf):
%   F(T1,T2) = forward + fwd_factor * f_{T1,T1}
%   S(T2)    = forward + fwd_factor * f_{T1,T1} + (f_{T2,T2} - f_{T1,T2})
% with fwd_factor = B(0,T1)/B(0,T2) = exp(int_{T1}^{T2} r_s ds).
%
% Per-model CDF engine:
%   * 'MA'      -> delegates to pricing_fwd_start_MA_MC (FA_simulation: IA base
%                  leg [0,T1] + FA increment leg [T1,T2]).
%   * 'AB','GL' -> samples both legs by inverting the Lewis-FFT digital CDF
%                  (lewis_FFT_digital) of the Lemma-2 rescaled increment CF
%                  (cf_increment_AB / cf_increment_GL), then PCHIP-inverts it
%                  (simulate_from_cdf). This mirrors the old price_fwd_start_MC
%                  but replaces ccdf_increment_FFT with lewis_FFT_digital.
%
% INPUTS
%   model    : 'MA' | 'AB' | 'GL'
%   params   : model parameter column vector
%                MA, GL -> [alpha; beta] ,  AB -> [k; eta]
%   sigma_T1 : FULL scale factor at the reset date T1, i.e.
%              (sigma_ATM(T1)/I_0) * sqrt(T1).  Pass 0 if T1 = 0.
%   sigma_T2 : FULL scale factor at maturity T2, i.e.
%              (sigma_ATM(T2)/I_0) * sqrt(T2).
%   N_sim    : number of MC paths
%   M, dz    : FFT grid parameters (FA_simulation for MA, Lewis-FFT for AB/GL)
%   forward  : F(0, T2)
%   B_0_t1   : B(0, T1)   (Lemma 2 forward rescaling)
%   B_0_t2   : B(0, T2)   (discounting; payoff received at T2)
%   strike   : K2 multiplier(s) (scalar or vector of strike multipliers)
%
% Reproducibility: this function does NOT seed the RNG; set rng(...) once in
% the caller so every MC stage shares a single deterministic stream.
%
% OUTPUTS
%   price : B(0,T2) * E[(S(T2) - K2*F(T1,T2))_+]   (1 x N_K row)
%   CI    : 95% confidence interval, [lower; upper]   (2 x N_K)
%   diag  : struct with the cleaned conditional CDF (x_cond, cdf_cond) and the
%           simulated increment W (for plotting / FFT cross-checks). Empty for
%           the MA branch.

    % --- Options Initialization ---
    if nargin < 12 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'verbose')
        opts.verbose = false;
    end
    if ~isfield(opts, 'plot')
        opts.plot = false;
    end
    % ------------------------------


    if N_sim == 0
        % Accuracy sizing decoupled into size_Nsim_MC and cached on disk
        sig_inputs = [params(:); sigma_T1; sigma_T2; forward; B_0_t1; B_0_t2; ...
                      strike(:); M; dz];
        N_sim = size_Nsim_MC( ...
            @(Np) nth_out(4, @pricing_fwd_start_MC, model, params, sigma_T1, sigma_T2, ...
                          Np, M, dz, forward, B_0_t1, B_0_t2, strike), ...
            sprintf('fwdstart_%s', upper(model)), sig_inputs, struct('ref', forward));
    end

    % Lemma 2
    fwd_factor = B_0_t1 / B_0_t2;

    switch upper(model)
        % -----------------------------------------------------------------
        case 'MA'   % Minimal Additive: dedicated FA_simulation engine
        % -----------------------------------------------------------------
            [price, CI, sigma] = pricing_fwd_start_MA_MC(forward, strike, B_0_t2, ...
                              N_sim, M, dz, [sigma_T1; sigma_T2], ...
                              params(1), params(2), fwd_factor, opts);
            diag = struct([]);

        % -----------------------------------------------------------------
        case {'AB', 'GL'}   % Additive Bachelier / Generalized Logistic
        % -----------------------------------------------------------------
            % Lemma-2 rescaled increment CF, signature (u, params, scale, fwd):
            %   phi_W(u) = phi_T2(u) / phi_T1(fwd_factor * u).
            if strcmpi(model, 'AB')
                cf_inc = @cf_increment_AB;
            else
                cf_inc = @cf_increment_GL;
            end
            % Wrap to the lewis_FFT_digital interface cf(u, params, scale).
            cf_h = @(u, p, sc) cf_inc(u, p, sc, fwd_factor);

            % --- marginal CDF (0 -> T1): scale(1)=0 collapses to phi_T1 -----
            [cdf_1, x_1] = lewis_FFT_digital(cf_h, M, dz, params, ...
                               [0; sigma_T1], true, upper(model), true, opts.plot);
            Z1 = simulate_from_cdf(cdf_1, x_1, true, N_sim);

            % --- conditional CDF (T1 -> T2) under Lemma 2 -------------------
            [cdf_2, x_2] = lewis_FFT_digital(cf_h, M, dz, params, ...
                               [sigma_T1; sigma_T2], true, upper(model), true, opts.plot);
            W = simulate_from_cdf(cdf_2, x_2, true, N_sim);

            % --- Lemma 2 reconstruction of forward / spot ------------------
            F_T1_T2 = forward + fwd_factor * Z1;       % F(T1, T2)
            S_T2    = forward + fwd_factor * Z1 + W;   % S(T2) = F(T2, T2)

            % --- payoff and MC estimator (vectorised over strikes) ---------
            K       = strike(:)';                      % 1 x N_K row
            payoff  = max(S_T2 - K .* F_T1_T2, 0);     % N_sim x N_K
            discounted_payoff  = B_0_t2 * payoff;
            [price, sigma, CI, ~] = normfit(discounted_payoff);

            diag = struct('x_cond', x_2, 'cdf_cond', cdf_2, 'W', W);

        otherwise
            error('pricing_fwd_start_MC:unknownModel', ...
                  'Unknown model "%s". Use ''MA'', ''AB'' or ''GL''.', model);
    end
end
