function [price, CI] = pricing_fwd_start_MC(forward, K, df, N_sim, M, dz, sigmat, alpha_MA, beta_MA)
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
%
%   OUTPUTS:
%       price     : The discounted expected payoff of the option
%       CI        : 2-element vector containing the [lower, upper] 95% Confidence Interval

    % =========================================================================
    % STEP 0: INPUT VALIDATION & INITIALIZATION
    % =========================================================================
    
    % =========================================================================
    % STEP 1: INTERNAL PARAMETER COMPUTATION
    % =========================================================================
    % 1.1 Asymmetry factor
    gamma_MA = (1/alpha_MA) - (1/beta_MA);
   
    % 1.4 Tail decay parameters for t1 (s) and t2 (t)
    % We strictly extract indices 2 and 3 to avoid division by zero at t=0.
    ps_plus  = beta_MA  / sigmat(1);
    ps_minus = alpha_MA / sigmat(1);
    
    pt_plus  = beta_MA  / sigmat(2);
    pt_minus = alpha_MA / sigmat(2);
    
    % 1.5 Drift vector calculation
    % drift(1): from 0 to t1. drift(2): from t1 to t2.
    drift_0_t1  = gamma_MA * (sigmat(1) - 0);
    drift_t1_t2 = gamma_MA * (sigmat(2) - sigmat(1));
    
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
    ft0t1 = FA_simulation(N_sim, M, dz, shift_pos_1, shift_neg_1, drift_0_t1, ...
                          ps_plus, ps_minus, 0, 0, 1, 'infinite', 1); 

    % =========================================================================
    % STEP 3: SIMULATE THE INCREMENT [t1 to t2]
    % =========================================================================
    % MA increment process -> 'finite' activity.
    % We pass both the target parameters (pt) and the subtracted parameters (ps).
    ft1t2 = FA_simulation(N_sim, M, dz, shift_pos_2, shift_neg_2, drift_t1_t2, ...
                          pt_plus, pt_minus, ps_plus, ps_minus, 1, 'finite', 1);

    % =========================================================================
    % STEP 4: ASSEMBLE PATHS & PAYOFF
    % =========================================================================
    % Ft1 represents the simulated underlying price at the forward start date (t1).
    Ft1 = forward + ft0t1; % Size: (N_sim, 1)
    
    % St2 represents the simulated final price at maturity (t2).
    St2 = Ft1 + ft1t2;     % Size: (N_sim, 1)
    
    % Ensure K is a row vector for implicit expansion
    K = K(:)'; % Size: (1, N_K)
    
    % Compute the payoff matrix. Size will be (N_sim, N_K).
    % St2 and Ft1 are expanded horizontally, K is expanded vertically.
    payoff = max(St2 - K .* Ft1, 0); 
    
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


