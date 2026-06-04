function [price, CI, sigma] = pricing_fwd_start_MA_MC(forward, K, df, N_sim, M, dz, sigmat, alpha_MA, beta_MA, fwd_factor, opts)
% PRICING_FWD_START_MA_MC Monte Carlo price of a forward-start option under MA.
%
%   This function encapsulates the Minimal Additive (MA) structural parameter
%   generation. From the core model parameters and the scale vector it computes
%   all internal quantities (drifts, tail decays, damping shifts) and simulates
%   the path in two legs: an Infinite-Activity base over [0, t1] plus a
%   Finite-Activity increment over [t1, t2] (see FA_simulation).
%
%   INPUTS:
%       forward   : forward F(0, t2)
%       K         : strike multiplier(s) (e.g. 1.0 for ATM forward start)
%       df        : discount factor B(0, t2) (payoff received at t2)
%       N_sim     : number of Monte Carlo paths
%       M, dz     : FFT grid parameters for the FA_simulation step
%       sigmat    : 2-element scale vector [sigma_t1*sqrt(t1); sigma_t2*sqrt(t2)]
%       alpha_MA  : alpha parameter of the MA process (governs left tail)
%       beta_MA   : beta parameter of the MA process (governs right tail)
%       fwd_factor: (optional, default 1) Lemma 2 rescaling B(0,t1)/B(0,t2)
%                   so that F(t1,t2) = forward + fwd_factor*f_{t1,t1}.
%       opts      : (optional) struct with fields .verbose and .plot
%
%   OUTPUTS:
%       price     : discounted expected payoff (1 x numel(K))
%       CI        : 95% confidence interval for the mean, [lower; upper]
%       sigma     : sample standard deviation of the discounted payoff

    % =========================================================================
    % STEP 0: INPUT VALIDATION & INITIALIZATION
    % =========================================================================
    % --- Options Initialization ---
    if nargin < 11 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'verbose')
        opts.verbose = false;
    end
    if ~isfield(opts, 'plot')
        opts.plot = false;
    end
    % ------------------------------

    % Lemma 2 
    if nargin < 10 || isempty(fwd_factor)
        fwd_factor = 1;
    end

    % if N_sim == 0
    %     % Run a pilot simulation to estimate standard deviation (e.g., 1000 paths)
    %     [~, ~, ~, sigma_est] = pricing_fwd_start_MA_MC(forward, K, df, 1000, M, dz, ...
    %         sigmat, alpha_MA, beta_MA, fwd_factor);
    %     target_error = 10 * 1e-4; 
    %     N_sim = min(ceil(((1.96 * sigma_est) / target_error)^2), 1e8);
    % 
    %     % Verification print
    %     fprintf('--- PILOT SIMULATION ---\n');
    %     fprintf('Estimated Std Dev: %.4f\n', sigma_est);
    %     fprintf('Target Error:      10 bps (%.4f)\n', target_error);
    %     fprintf('Required N_sim:    %d\n', N_sim);
    %     fprintf('------------------------\n');
    % end
    

    % =========================================================================
    % STEP 1: INTERNAL PARAMETER COMPUTATION
    % =========================================================================
    % 1.1 Asymmetry factor
    gamma_MA = (1/alpha_MA) - (1/beta_MA);

    % Model parameters as a column vector 
    params = [alpha_MA; beta_MA];
   
    % 1.4 Tail decay parameters for t1 (s) and t2 (t).
    sigmat_s = fwd_factor * sigmat(1);

    ps_plus  = beta_MA  / sigmat_s;
    ps_minus = alpha_MA / sigmat_s;

    pt_plus  = beta_MA  / sigmat(2);
    pt_minus = alpha_MA / sigmat(2);

    % 1.5 Drift vector calculation
    drift_0_t1  = gamma_MA * (sigmat(1) - 0);
    drift_t1_t2 = gamma_MA * (sigmat(2) - sigmat_s);
    
    % 1.6 FFT Damping shifts (using half of the poles bounds)
    shift_pos_1 =  ps_minus / 2;
    shift_neg_1 = -ps_plus  / 2;
    
    shift_pos_2 =  pt_minus / 2;
    shift_neg_2 = -pt_plus  / 2;
    if opts.verbose
        fprintf('\n==================================================\n');
        fprintf('   FA_SIMULATION RUN DIAGNOSTICS                  \n');
        fprintf('==================================================\n');
        fprintf('General Grid Parameters:\n');
        fprintf('  -> N_sim       : %d\n', N_sim);
        fprintf('  -> M (Grid Size): %d (N = %d points)\n', M, 2^M);
        fprintf('  -> dz          : %.6f\n', dz);
        
        fprintf('\nFirst Leg [0 to t1] - INFINITE ACTIVITY:\n');
        fprintf('  -> drift_0_t1  : %.6f\n', drift_0_t1);
        fprintf('  -> shift_pos_1 : %.6f\n', shift_pos_1);
        fprintf('  -> shift_neg_1 : %.6f\n', shift_neg_1);
        fprintf('  -> ps_plus     : %.6f\n', ps_plus);
        fprintf('  -> ps_minus    : %.6f\n', ps_minus);
        
        fprintf('\nSecond Leg [t1 to t2] - FINITE ACTIVITY:\n');
        fprintf('  -> drift_t1_t2 : %.6f\n', drift_t1_t2);
        fprintf('  -> shift_pos_2 : %.6f\n', shift_pos_2);
        fprintf('  -> shift_neg_2 : %.6f\n', shift_neg_2);
        fprintf('  -> pt_plus     : %.6f\n', pt_plus);
        fprintf('  -> pt_minus    : %.6f\n', pt_minus);
        fprintf('  -> ps_plus (sub): %.6f\n', ps_plus);
        fprintf('  -> ps_minus(sub): %.6f\n', ps_minus);
        fprintf('==================================================\n\n');
    end
    % =========================================================================
    % STEP 2: SIMULATE THE FIRST LEG [0 to t1]
    % =========================================================================

    ft0t1 = FA_simulation(N_sim, M, dz, drift_0_t1, ...
                          ps_plus, ps_minus, 0, 0, 1, 'infinite', 1, ...
                          params, sigmat(1), opts.plot);

    % =========================================================================
    % STEP 3: SIMULATE THE INCREMENT [t1 to t2]
    % =========================================================================

    ft1t2 = FA_simulation(N_sim, M, dz, drift_t1_t2, ...
                          pt_plus, pt_minus, ps_plus, ps_minus, 1, 'finite', 1, ...
                          params, [sigmat_s; sigmat(2)], opts.plot);

    % =========================================================================
    % STEP 4: ASSEMBLE PATHS & PAYOFF
    % =========================================================================

    F_T1_T2 = forward + fwd_factor * ft0t1; % Size: (N_sim, 1)

    St2 = forward + fwd_factor * ft0t1 + ft1t2;     % Size: (N_sim, 1)

    K = K(:)'; % Size: (1, N_K)

    payoff = max(St2 - K .* F_T1_T2, 0);
    discounted_payoff = df * payoff;
    
    % =========================================================================
    % STEP 5: FINAL PRICING AND CONFIDENCE INTERVAL (95%)
    % =========================================================================
    [price, sigma, CI, ~] = normfit(discounted_payoff);

end


