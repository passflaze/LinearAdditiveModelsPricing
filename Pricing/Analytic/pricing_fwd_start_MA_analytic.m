function price = pricing_fwd_start_MA_analytic(alpha_MA, beta_MA, sigmat, df, K2, F_t0_t2, fwd_factor)
% PRICING_FWD_START_ANALYTIC Computes the analytical price of a forward start
% option under the Minimal Additive model. Supports vectorized strikes (K2).
%
% Inputs:
%   alpha_MA, beta_MA : Tail parameters
%   sigmat            : Vector of integrated volatilities [sigma*sqrt(t1), sigma*sqrt(t2)]
%   df                : Total discount factor from t0 to t2 ( B(t0,t2) )
%   K2                : Proportional strike multiplier (Scalar or Vector)
%   F_t0_t2           : Forward price at t0 expiring at t2
%   fwd_factor        : (optional, default 1) Lemma 2 (Forward.pdf) rescaling
%                       B(0,T1)/B(0,T2). With the change of variable
%                       X' = fwd_factor*f_{T1,T1} = f_{T1,T2}, the forward-start
%                       payoff and the joint law map EXACTLY onto the
%                       fwd_factor = 1 formulation with the s leg scale rescaled
%                       to fwd_factor*sigmat(1) (MA marginals scale linearly).
%                       Hence we simply rescale sigmat(1) below.

    if nargin < 7 || isempty(fwd_factor)
        fwd_factor = 1;
    end
    % Lemma 2 rescaling of the subtracted (s = t1) leg. sigmat(2) is unchanged.
    sigmat(1) = fwd_factor * sigmat(1);

    % =========================================================================
    % STEP 1: COMPUTE STRUCTURAL PARAMETERS
    % =========================================================================
    gamma_MA = (1/alpha_MA) - (1/beta_MA);
    ps_plus  = beta_MA  / sigmat(1); % T1 parameters (s)
    ps_minus = alpha_MA / sigmat(1);
    pt_plus  = beta_MA  / sigmat(2); % T2 parameters (t)
    pt_minus = alpha_MA / sigmat(2);
    
    % Drift (Delta_mu) and Marginal Shift (mu_T1)
    drift_t1_t2 = gamma_MA * (sigmat(2) - sigmat(1)); % Delta mu
    drift_t0_t1 = gamma_MA * sigmat(1);               % mu_T1
    
    % =========================================================================
    % STEP 2: COMPUTE ANALYTICAL CONSTANTS
    % =========================================================================
    c = (pt_plus * pt_minus) / (ps_plus * ps_minus);
    denom = pt_plus + pt_minus;
    
    A_minus = c * ((ps_minus - pt_minus) * (ps_plus + pt_minus)) / denom;
    A_plus  = c * ((ps_plus - pt_plus) * (ps_minus + pt_plus)) / denom;
    
    % Forward Start Specific Constants (Vectors)
    K_star = K2 - 1;
    
    % Inizializzazione degli array di output e appoggio
    price_undiscounted = zeros(size(K2));
    x_hat = zeros(size(K2));
    
    % --- EDGE CASE HANDLING (K2 = 1) ---
    idx_atm = (K_star == 0);
    if any(idx_atm)
        if drift_t1_t2 > 0 % Strike 0 < Delta_mu (Inner ITM)
            atm_val = (A_minus / pt_minus^2) * exp(-pt_minus * drift_t1_t2);
        else               % Strike 0 >= Delta_mu (Inner OTM)
            atm_val = (A_plus / pt_plus^2) * exp(pt_plus * drift_t1_t2);
        end
        price_undiscounted(idx_atm) = atm_val;
    end
    
    % --- NON-ATM HANDLING ---
    idx_non_atm = (K_star ~= 0);
    mu_t1 = drift_t0_t1;
    
    % Calcoliamo x_hat solo dove K_star non è zero per evitare divisioni per zero
    x_hat(idx_non_atm) = (drift_t1_t2 ./ K_star(idx_non_atm)) - F_t0_t2;
    
    % Substitution Constants (Vettorizzate)
    C_lin = K_star .* F_t0_t2;
    C_in  = (A_minus / pt_minus^2) .* exp(pt_minus .* (K_star .* F_t0_t2 - drift_t1_t2));
    C_out = (A_plus  / pt_plus^2)  .* exp(-pt_plus .* (K_star .* F_t0_t2 - drift_t1_t2));
    
    % Density Scale Constants (Scalari)
    D_L = (ps_plus * ps_minus) / (ps_plus + ps_minus) * exp(-ps_minus * mu_t1);
    D_R = (ps_plus * ps_minus) / (ps_plus + ps_minus) * exp(ps_plus * mu_t1);
    
    % =========================================================================
    % STEP 3: PRIMITIVE FUNCTIONS EVALUATION (Psi_j)
    % =========================================================================
    % Aggiornate per accettare vettori (x, k, clin, cin, cout) e operare element-wise
    
    Psi1 = @(x, k, clin, cin) D_L .* ( - (clin ./ ps_minus) .* exp(ps_minus .* x) ...
                        - k .* ((ps_minus .* x - 1) ./ ps_minus^2) .* exp(ps_minus .* x) ...
                        + (cin ./ (ps_minus + pt_minus .* k)) .* exp((ps_minus + pt_minus .* k) .* x) );
                    
    Psi2 = @(x, k, clin, cin) D_R .* (   (clin ./ ps_plus) .* exp(-ps_plus .* x) ...
                        + k .* ((ps_plus .* x + 1) ./ ps_plus^2) .* exp(-ps_plus .* x) ...
                        - (cin ./ (ps_plus - pt_minus .* k)) .* exp(-(ps_plus - pt_minus .* k) .* x) );
                    
    Psi3 = @(x, k, cout) D_L .* (   (cout ./ (ps_minus - pt_plus .* k)) .* exp((ps_minus - pt_plus .* k) .* x) );
    
    Psi4 = @(x, k, cout) -D_R .* (  (cout ./ (ps_plus + pt_plus .* k)) .* exp(-(ps_plus + pt_plus .* k) .* x) );
    
    % =========================================================================
    % STEP 4: PIECEWISE OPTION PRICING (LOGICAL MASKS)
    % =========================================================================
    
    % Creazione delle maschere per i 4 regimi
    idx_1A = (K_star > 0) & (x_hat < mu_t1) & idx_non_atm;
    idx_1B = (K_star > 0) & (x_hat >= mu_t1) & idx_non_atm;
    idx_2A = (K_star < 0) & (x_hat < mu_t1) & idx_non_atm;
    idx_2B = (K_star < 0) & (x_hat >= mu_t1) & idx_non_atm;
    
    % Valutazione Case 1A
    if any(idx_1A)
        price_undiscounted(idx_1A) = Psi1(x_hat(idx_1A), K_star(idx_1A), C_lin(idx_1A), C_in(idx_1A)) ...
                                   + Psi3(mu_t1, K_star(idx_1A), C_out(idx_1A)) ...
                                   - Psi3(x_hat(idx_1A), K_star(idx_1A), C_out(idx_1A)) ...
                                   - Psi4(mu_t1, K_star(idx_1A), C_out(idx_1A));
    end
    
    % Valutazione Case 1B
    if any(idx_1B)
        price_undiscounted(idx_1B) = Psi1(mu_t1, K_star(idx_1B), C_lin(idx_1B), C_in(idx_1B)) ...
                                   + Psi2(x_hat(idx_1B), K_star(idx_1B), C_lin(idx_1B), C_in(idx_1B)) ...
                                   - Psi2(mu_t1, K_star(idx_1B), C_lin(idx_1B), C_in(idx_1B)) ...
                                   - Psi4(x_hat(idx_1B), K_star(idx_1B), C_out(idx_1B));
    end
    
    % Valutazione Case 2A
    if any(idx_2A)
        price_undiscounted(idx_2A) = Psi3(x_hat(idx_2A), K_star(idx_2A), C_out(idx_2A)) ...
                                   + Psi1(mu_t1, K_star(idx_2A), C_lin(idx_2A), C_in(idx_2A)) ...
                                   - Psi1(x_hat(idx_2A), K_star(idx_2A), C_lin(idx_2A), C_in(idx_2A)) ...
                                   - Psi2(mu_t1, K_star(idx_2A), C_lin(idx_2A), C_in(idx_2A));
    end
    
    % Valutazione Case 2B
    if any(idx_2B)
        price_undiscounted(idx_2B) = Psi3(mu_t1, K_star(idx_2B), C_out(idx_2B)) ...
                                   + Psi4(x_hat(idx_2B), K_star(idx_2B), C_out(idx_2B)) ...
                                   - Psi4(mu_t1, K_star(idx_2B), C_out(idx_2B)) ...
                                   - Psi2(x_hat(idx_2B), K_star(idx_2B), C_lin(idx_2B), C_in(idx_2B));
    end
    
    % Applichiamo il fattore di sconto totale (B(t0, t2)) vettorizzato
    price = df .* price_undiscounted;
end