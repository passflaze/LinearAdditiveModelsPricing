% =========================================================================
%  HEDGING MODULE — SHARED DATA CONTRACTS  (Ex.6 Risk Management, AB model)
% =========================================================================
%  This file is documentation only (not executed). It fixes the struct
%  interfaces shared between the two work streams so they can be developed
%  in parallel. Do NOT change a field name without telling the other owner.
%
%  Conventions (project-wide): params_AB = [k; eta] (column);
%  scale_factor = [scale_t1, scale_t2] = sigma_t .* sqrt(yf), with
%  sigma_t = sigma_ATM / I0_AB(0, params_AB). Vanillas/future on T2.
%
%  -----------------------------------------------------------------------
%  mkt    : market snapshot for one exotic
%    .forward   F(t0,T2)              (scalar)
%    .K1        compound strike       (scalar, CoC/PoP only)
%    .K2        inner strike = F(T2,T2)
%    .df        [B(t0,T1), B(t0,T2)]  (1x2)
%    .Kc, .Kp   strikes of the hedge call / put (fixed from t0).
%               MUST be at different moneyness (Kc != Kp): a same-strike
%               call/put pair is collinear in (gamma, vega) -> build_hedge_AB
%               singular. run_ex6 snaps an OTM call and a further-OTM put.
%
%  contract : FIXED terms of one exotic+hedge, carried across snapshots
%    .K1 .K2 .Kc .Kp   strikes (never change)
%    .E1 .E2           calendar expiries T1, T2 (datetime), matched by DATE
%                      inside recalibrate_AB (indices shift as maturities roll)
%
%  mc     : Monte Carlo settings (exotic pricers)
%    .N_sim .M .dz .N_grid .seed      (.seed -> common random numbers)
%
%  bumps  : finite-difference steps (SAME for exotic and vanillas)
%    .dF    absolute forward bump ($)              e.g. 0.5   (delta, gamma)
%    .dSig  relative ATM-vol bump                  e.g. 1e-2  (scale-bump vega)
%    .dvol  absolute sticky-strike vol shift       e.g. 0.5   (recalib vega)
%
%  VEGA MODES (run_ex6 vega_mode):
%    'scale'   analytic bump of the normalized scale: self-similar /
%              sticky-moneyness move, (k,eta) FIXED. No re-calibration, exact,
%              low noise. Units: per RELATIVE vol shift (dSig).
%    'recalib' sticky-strike PARALLEL shift of the implied-vol surface + AB
%              re-calibration (shock_recalibrate_AB): (k,eta) MOVE,
%              model-consistent. Units: per ABSOLUTE Bachelier vol (dvol).
%    The two answer different smile-dynamics questions; run_ex6 reports both
%    plus the (k,eta) drift, and hedges with the recalib (sticky-strike) vega.
%
%  greeks : output of greeks_exotic_AB / greeks_vanilla_AB
%    .price .delta .gamma .vega
%
%  positions : output of build_hedge_AB (per 1 unit of exotic)
%    .nC .nP .nF
%
%  costRule : bid-ask rule (cost on the notional F)
%    .fut_bp = 1   .opt_bp = 4
%
%  state  : full per-date snapshot used by hedge_backtest
%    .params_AB .scale_factor .mkt
%
%  P&L sign convention (hedge_backtest): book = LONG 1 exotic, SHORT the
%  replicating hedge. The future enters the P&L as MARK-TO-MARKET,
%  nF*(F_t - F_{t-1}), NOT as nF*F (a future has zero entry cost).
%
%  -----------------------------------------------------------------------
%  FILE MAP
%    price_exotic_AB        dispatch CoC/PoP/Chooser MC (+ CRN)
%    vanilla_AB_price       single AB vanilla (call/put), shared building block
%    greeks_exotic_AB       delta/gamma (forward bump) + scale-bump vega
%    greeks_vanilla_AB      delta/gamma/vega of call/put/future
%    shock_recalibrate_AB   sticky-strike IV shock + AB re-calibration
%    vega_recalib_sticky    model-consistent vega (+/- shocks shared, CRN)
%    build_hedge_AB         delta-gamma-vega hedge solve (+ cond guard)
%    hedging_cost           bid-ask transaction cost
%    hedge_backtest         next-Tuesdays P&L / cost backtest
%    run_ex6                orchestrator (vega_mode = scale|recalib|both)
%
%  OWNERSHIP
%    Persona A : price_exotic_AB, vanilla_AB_price, greeks_exotic_AB,
%                greeks_vanilla_AB, shock_recalibrate_AB, vega_recalib_sticky
%    Persona B : build_hedge_AB, hedging_cost, hedge_backtest, run_ex6
%  Sync points : the structs above.
% =========================================================================
