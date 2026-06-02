function out = vega_recalib_sticky(exotics, market, params0, mkt, mc, iT1, iT2, dvol)
% VEGA_RECALIB_STICKY  Model-consistent vega via sticky-strike re-calibration.
%
%   Central difference of the exotic and hedge-vanilla prices under a +/- dvol
%   sticky-strike shock of the implied-vol surface, with the AB model
%   re-calibrated at each shock (shock_recalibrate_AB). The two heavy
%   re-calibrations (+dvol, -dvol) are done ONCE and shared across all exotics
%   and both vanillas. CRN (mc.seed reset inside price_exotic_AB) cancels the
%   MC noise in the exotic differences.
%
%   Units: vega is per ABSOLUTE Bachelier-vol unit (d/dvol). The hedge ratios
%   are consistent because the exotic AND the vanillas are differenced with
%   the SAME dvol. Only the vega component differs from the scale-bump greeks;
%   delta and gamma (pure forward bumps) are unchanged and taken from there.
%
% INPUTS:
%   exotics - (cell) subset of {'CoC','PoP','Chooser'}
%   market  - (struct) run_ex2 output (see shock_recalibrate_AB)
%   params0 - (2x1) base AB params [k; eta]
%   mkt     - (struct) exotic/hedge market snapshot (forward, K1, K2, Kc, Kp, df)
%   mc      - (struct) MC settings (price_exotic_AB)
%   iT1,iT2 - (scalar) maturity indices of the (T1,T2) window
%   dvol    - (scalar) sticky-strike vol bump (absolute Bachelier vol)
%
% OUTPUT:
%   out - (struct)
%         .exotic.(name)  recalib vega of each exotic
%         .call .put      recalib vega of the two hedge vanillas
%         .drift          [k_up k_dn; eta_up eta_dn] - params0 (shape response)
%         .params_up .params_dn  recalibrated params at +/- dvol
%         .sig_up .sig_dn        shocked ATM vols (at iT2) for reporting

    % --- Two shared re-calibrations --------------------------------------
    s_up = shock_recalibrate_AB(market, params0, +dvol);
    s_dn = shock_recalibrate_AB(market, params0, -dvol);

    scale_up = scale_at(s_up, market.yf, iT1, iT2);
    scale_dn = scale_at(s_dn, market.yf, iT1, iT2);

    % --- Vanilla vegas (call at Kc, put at Kp) ---------------------------
    cC_up = vanilla_AB_price(mkt.forward, mkt.Kc, scale_up(2), mkt.df(2), s_up.params_AB, true);
    cC_dn = vanilla_AB_price(mkt.forward, mkt.Kc, scale_dn(2), mkt.df(2), s_dn.params_AB, true);
    cP_up = vanilla_AB_price(mkt.forward, mkt.Kp, scale_up(2), mkt.df(2), s_up.params_AB, false);
    cP_dn = vanilla_AB_price(mkt.forward, mkt.Kp, scale_dn(2), mkt.df(2), s_dn.params_AB, false);

    out.call = (cC_up - cC_dn) / (2*dvol);
    out.put  = (cP_up - cP_dn) / (2*dvol);

    % --- Exotic vegas (shared shocks, CRN) -------------------------------
    out.exotic = struct();
    for ee = 1:numel(exotics)
        nm   = exotics{ee};
        V_up = price_exotic_AB(nm, s_up.params_AB, scale_up, mkt, mc);
        V_dn = price_exotic_AB(nm, s_dn.params_AB, scale_dn, mkt, mc);
        out.exotic.(nm) = (V_up - V_dn) / (2*dvol);
    end

    % --- Diagnostics ------------------------------------------------------
    out.params_up = s_up.params_AB;
    out.params_dn = s_dn.params_AB;
    out.drift     = [s_up.params_AB, s_dn.params_AB] - params0(:);
    out.sig_up    = s_up.sigma_ATM(iT2);
    out.sig_dn    = s_dn.sigma_ATM(iT2);
end

% =========================================================================
function scale = scale_at(state, yf, iT1, iT2)
% AB normalized scale [scale_t1, scale_t2] from a shocked state.
    I0  = I0_AB(0, state.params_AB);
    sig = state.sigma_ATM / I0;
    scale = [sig(iT1)*sqrt(yf(iT1)), sig(iT2)*sqrt(yf(iT2))];
end
