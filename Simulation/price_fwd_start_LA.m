function [price, CI, ref_atm] = price_fwd_start_LA(spec, T1, T2, ...
                                       sigma_T1, sigma_T2, forward, ...
                                       B_0_T1, B_0_T2, K2, N_sim)
%PRICE_FWD_START_LA  Monte Carlo price of the forward-start option
%   payoff(T2) = [ S(T2) - K2 * F(T1,T2) ]_+    (received at T2)
%   under a Linear Additive model (AB or GL), with the analytical/FFT ATM
%   reference. Unified pricer for the two "twin" models: it mirrors the
%   structure of price_AB_MC but drives the model through the SPEC struct and
%   is vectorized over the strike multiplier K2.
%
%   Strategy (additivity of the AB/GL increments + Lemma 2, Forward.pdf):
%     1. Z1 ~ f_{T1,T1}            from the marginal CDF (0 -> T1)
%     2. W  ~ f_{T2,T2}-f_{T1,T2}  from the conditional CDF (T1 -> T2)
%     3. Lemma 2: f_{T1,T2} = fwd_factor * f_{T1,T1}, fwd_factor = B(0,T1)/B(0,T2)
%          F(T1,T2) = forward + fwd_factor * Z1
%          S(T2)    = forward + fwd_factor * Z1 + W
%
%   Note: for K2 = 1 the payoff collapses to max(W,0) (independent of Z1), so the
%   exact reference is E[max(W,0)] = int_0^inf (1 - F_W(x)) dx, returned in ref_atm.
%
%   INPUTS
%     spec      : model spec struct (la_model_spec)
%     T1, T2    : reset and maturity year fractions
%     sigma_T1  : model scale at T1 (= sigma_ATM(T1)/I0)
%     sigma_T2  : model scale at T2 (= sigma_ATM(T2)/I0)
%     forward   : F(0,T2) initial forward at T2
%     B_0_T1    : B(0,T1) discount factor (Lemma-2 forward rescaling)
%     B_0_T2    : B(0,T2) discount factor (payoff received at T2)
%     K2        : strike multiplier(s), scalar or vector
%     N_sim     : number of MC paths
%
%   OUTPUTS
%     price   : 1 x nK discounted MC prices
%     CI      : 2 x nK 95% confidence interval [lower; upper]
%     ref_atm : scalar exact ATM (K2 = 1) reference price (survival integral)

    fwd_factor = B_0_T1 / B_0_T2;          % = exp(int_{T1}^{T2} r_s ds)

    % --- (1) marginal increment 0 -> T1 --------------------------------------
    Z1 = simulate_LA_increment(spec, 0, T1, 0, sigma_T1, 1, N_sim);

    % --- (2) conditional increment T1 -> T2 (keep CDF for the ATM reference) --
    [W, x_fine, cdf_fine] = simulate_LA_increment(spec, T1, T2, ...
                                       sigma_T1, sigma_T2, fwd_factor, N_sim);

    % --- (3) reconstruct forward / spot under Lemma 2 ------------------------
    F1 = forward + fwd_factor * Z1;        % F(T1,T2)   (N_sim x 1)
    S2 = F1 + W;                           % S(T2)      (N_sim x 1)

    % --- payoff vectorized over strikes -------------------------------------
    K = K2(:).';                           % 1 x nK
    payoff = max(S2 - K .* F1, 0);         % N_sim x nK

    price     = B_0_T2 * mean(payoff, 1);          % 1 x nK
    std_err   = B_0_T2 * std(payoff, 0, 1) / sqrt(N_sim);
    CI        = [price - 1.96*std_err; price + 1.96*std_err];

    % --- exact ATM (K2 = 1) reference: E[max(W,0)] = int_0^inf (1-F_W) dx -----
    % x_fine may not hit 0 exactly: interpolate F_W(0) so the integral starts at 0.
    F0       = interp1(x_fine, cdf_fine, 0, 'pchip');
    mask_pos = x_fine > 0;
    ref_atm  = B_0_T2 * trapz([0; x_fine(mask_pos)], 1 - [F0; cdf_fine(mask_pos)]);
end
