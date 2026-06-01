function [price, CI] = pricing_fwd_start_MA_MC(forward, K, df, N_sim, M, dz, sigmat, alpha_MA, beta_MA, fwd_factor)
% PRICING_FWD_START_MC Computes the Monte Carlo price of a Forward Start Option.
%
%   This function completely encapsulates the Minimal Additive (MA) structural
%   parameter generation. It takes the core model parameters and a time vector,
%   computes all internal variables (drifts, tail decays, shifts), and simulates 
%   the paths in two steps (Infinite Activity base + Finite Activity increment).
%
%   INPUTS:
%       forward   : Initial forward/spot price at t=0
%       K         : Strike multiplier (e.g., 1.0 for ATM forward start)
%       df        : Discount factor from maturity (t2) to present (0)
%       N_sim     : Number of Monte Carlo paths
%       M, dz     : FFT grid parameters for the FA_simulation step
%       t_vec     : 3-element time vector [0, t1, t2] 
%       sigma_ATM : At-The-Money implied volatility
%       alpha_MA  : Alpha parameter of the MA process (governs left tail)
%       beta_MA   : Beta parameter of the MA process (governs right tail)
%       fwd_factor: (optional, default 1) Lemma 2 rescaling B(0,t1)/B(0,t2)
%                   so that F(t1,t2) = forward + fwd_factor*f_{t1,t1}.
%
%   OUTPUTS:
%       price     : The discounted expected payoff of the option
%       CI        : 2-element vector containing the [lower, upper] 95% Confidence Interval

    % =========================================================================
    % STEP 0: INPUT VALIDATION & INITIALIZATION
    % =========================================================================
    % Lemma 2 (Forward.pdf): f_{T1,T2} = fwd_factor * f_{T1,T1} with
    % fwd_factor = B(0,T1)/B(0,T2) = exp(int_{T1}^{T2} r_s ds). Default 1
    % recovers the zero-inter-reset-rate case (backward compatible).
    if nargin < 10 || isempty(fwd_factor)
        fwd_factor = 1;
    end

    % =========================================================================
    % STEP 1: INTERNAL PARAMETER COMPUTATION
    % =========================================================================
    % 1.1 Asymmetry factor
    gamma_MA = (1/alpha_MA) - (1/beta_MA);

    % Model parameters as a column vector [alpha; beta] for the unified engine
    % (lewis_FFT_digital / conditional_cf_MA_FA / cf_MA_IA).
    params = [alpha_MA; beta_MA];
   
    % 1.4 Tail decay parameters for t1 (s) and t2 (t).
    % Lemma 2 (Forward.pdf): the t1->t2 increment belongs to the SAME-expiry
    % forward F(.,T2), whose reset value is f_{T1,T2} = fwd_factor*f_{T1,T1}.
    % Hence the SUBTRACTED (s = t1) leg of the increment carries the rescaled
    % scale sigma_s = fwd_factor*sigmat(1), i.e. poles ps_+- = a,b/sigma_s (the
    % original poles divided by fwd_factor). This mirrors phi_t1(fwd_factor*u)
    % in cf_increment_AB/GL. The 0->t1 base leg below keeps the TRUE marginal
    % scale sigmat(1), since it samples f_{T1,T1} directly.
    sigmat_s = fwd_factor * sigmat(1);

    ps_plus  = beta_MA  / sigmat_s;
    ps_minus = alpha_MA / sigmat_s;

    pt_plus  = beta_MA  / sigmat(2);
    pt_minus = alpha_MA / sigmat(2);

    % 1.5 Drift vector calculation
    % drift(1): 0 -> t1 base leg, TRUE marginal scale sigmat(1).
    % drift(2): t1 -> t2 increment, Lemma-2 rescaled s leg (sigma_s).
    drift_0_t1  = gamma_MA * (sigmat(1) - 0);
    drift_t1_t2 = gamma_MA * (sigmat(2) - sigmat_s);
    
    % 1.6 FFT Damping shifts (using half of the poles bounds)
    shift_pos_1 =  ps_minus / 2;
    shift_neg_1 = -ps_plus  / 2;
    
    shift_pos_2 =  pt_minus / 2;
    shift_neg_2 = -pt_plus  / 2;

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
    % =========================================================================
    % STEP 2: SIMULATE THE FIRST LEG [0 to t1]
    % =========================================================================
    % Base process evaluated at time t1 -> 'infinite' activity.
    % Subtracted tail parameters are 0 because it starts from 0.
    % Unified engine: 'infinite' leg uses the SCALAR scale sigmat(1) (=t1),
    % since cf_MA_IA collapses scale_factors internally. delta_mu is ignored
    % here (the drift is already baked into cf_MA_IA).
    ft0t1 = FA_simulation(N_sim, M, dz, drift_0_t1, ...
                          ps_plus, ps_minus, 0, 0, 1, 'infinite', 1, ...
                          params, sigmat(1), []);

    % =========================================================================
    % STEP 3: SIMULATE THE INCREMENT [t1 to t2]
    % =========================================================================
    % MA increment process -> 'finite' activity.
    % We pass both the target parameters (pt) and the subtracted parameters (ps).
    % Unified engine: 'finite' increment leg uses the 2-VECTOR scale
    % [sigmat(1); sigmat(2)] so conditional_cf_MA_FA gets both s=t1 and t=t2
    % tails. delta_mu = drift_t1_t2 IS added back in this branch.
    ft1t2 = FA_simulation(N_sim, M, dz, drift_t1_t2, ...
                          pt_plus, pt_minus, ps_plus, ps_minus, 1, 'finite', 1, ...
                          params, [sigmat_s; sigmat(2)], []);

    % =========================================================================
    % STEP 4: ASSEMBLE PATHS & PAYOFF
    % =========================================================================
    % Strike reference: forward F(T1,T2) seen at the reset date t1. By Lemma 2
    % f_{T1,T2} = fwd_factor * f_{T1,T1}, hence F(T1,T2) = forward + fwd_factor*ft0t1.
    F_T1_T2 = forward + fwd_factor * ft0t1; % Size: (N_sim, 1)

    % St2 = forward + f_{T2,T2}, with f_{T2,T2} = f_{T1,T2} + (f_{T2,T2}-f_{T1,T2})
    %                                   = fwd_factor*f_{T1,T1} + W   (Lemma 2).
    % ft0t1 therefore enters with coefficient fwd_factor (the reset instantaneous
    % forward is f_{T1,T2}, NOT f_{T1,T1}), matching the AB/GL branch
    % S_T2 = forward + fwd_factor*Z1 + W in pricing_fwd_start_MC.
    St2 = forward + fwd_factor * ft0t1 + ft1t2;     % Size: (N_sim, 1)

    % Ensure K is a row vector for implicit expansion
    K = K(:)'; % Size: (1, N_K)

    % Compute the payoff matrix. Size will be (N_sim, N_K).
    % St2 and F_T1_T2 are expanded horizontally, K is expanded vertically.
    payoff = max(St2 - K .* F_T1_T2, 0);
    
    % =========================================================================
    % STEP 5: FINAL PRICING AND CONFIDENCE INTERVAL (95%)
    % =========================================================================
    % Calculate the mean and standard deviation along the first dimension (paths)
    payoff_mean = mean(payoff, 1); % Size: (1, N_K)
    payoff_std  = std(payoff, 0, 1); % Size: (1, N_K)
    
    % The final price is the discounted expected payoff
    price = df * payoff_mean;
    
    % Standard Error = Standard Deviation / sqrt(N)
    std_error = (df * payoff_std) / sqrt(N_sim);
    
    % 95% Confidence Interval bounds (Z-score = 1.96)
    CI_lower = price - 1.96 * std_error;
    CI_upper = price + 1.96 * std_error;
    
    CI = [CI_lower; CI_upper];

end


