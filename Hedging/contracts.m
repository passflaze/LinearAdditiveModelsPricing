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
%    .dF    absolute forward bump ($)         e.g. 0.5
%    .dSig  relative ATM-vol bump             e.g. 1e-2
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
%  OWNERSHIP
%    Persona A : price_exotic_AB, greeks_exotic_AB, greeks_vanilla_AB
%    Persona B : build_hedge_AB, hedging_cost, hedge_backtest, run_ex6
%  Sync points : the structs above.
% =========================================================================
