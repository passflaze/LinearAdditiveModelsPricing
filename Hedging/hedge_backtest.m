function res = hedge_backtest(type, contract, state0, positions0, tuesdays, dataOpts, mc, bumps, costRule)
% HEDGE_BACKTEST  Self-financing-style hedge check over the next Tuesdays.
%
%   The book is LONG 1 exotic and SHORT the replicating hedge produced by
%   build_hedge_AB (positions = quantities that replicate the exotic greeks;
%   we hold the NEGATIVE of them to neutralise a long exotic). On each
%   Tuesday the market is re-observed, the AB model is RE-CALIBRATED, the
%   exotic and the hedge instruments are repriced, the realised step P&L is
%   booked and the hedge is rebalanced (paying the bid-ask cost).
%
%   P&L convention (per step t_{k-1} -> t_k, hedge held = positions at t_{k-1}):
%       dV_exotic   = V_ex(t_k) - V_ex(t_{k-1})
%       dHedge      = nC*(C_k - C_{k-1}) + nP*(P_k - P_{k-1})
%                                        + nF*(F_k - F_{k-1})   <- future is MtM
%       PnL_unhedged = dV_exotic
%       PnL_hedged   = dV_exotic - dHedge        (long exotic, short hedge)
%   A well-hedged book leaves only theta + un-hedged skew/convexity in
%   PnL_hedged, while PnL_unhedged carries the full directional move.
%
% INPUTS:
%   type       - (string) 'CoC' | 'PoP' | 'Chooser'
%   contract   - (struct) FIXED contract terms (strikes never change):
%                  .K1 .K2 .Kc .Kp        option strikes
%                  .E1 .E2                calendar expiries T1, T2 (datetime)
%   state0     - (struct) inception state at t0:
%                  .params_AB .scale_factor .mkt   (see price_exotic_AB)
%   positions0 - (struct) .nC .nP .nF  hedge set at t0 (build_hedge_AB)
%   tuesdays   - (Nx1 datetime) value dates to test (e.g. next two Tuesdays)
%   dataOpts   - (struct) .callpath .putpath .expiryFile  (for run_ex2)
%   mc, bumps, costRule - see greeks_exotic_AB / hedging_cost
%
% OUTPUT:
%   res - (table) one row per date:
%         Date V_exotic PnL_unhedged PnL_hedged cost_step cost_cum
%         res_delta res_gamma res_vega

    % --- Inception marks (t0) --------------------------------------------
    prev = mark_state(type, state0, mc, bumps);
    prev_pos  = positions0;
    cost_cum  = hedging_cost(positions0, prev.F, costRule);   % set-up cost

    PnL_unh_cum = 0;
    PnL_hed_cum = 0;

    nT   = numel(tuesdays);
    rows = cell(nT, 1);

    for kk = 1:nT
        snap = tuesdays(kk);

        % --- Re-observe market & re-calibrate AB at this Tuesday ----------
        state_k = recalibrate_AB(snap, dataOpts, contract);
        cur     = mark_state(type, state_k, mc, bumps);

        % --- Realised step P&L (hedge held over the step = prev_pos) ------
        dV   = cur.Vex - prev.Vex;
        dHed = prev_pos.nC*(cur.C - prev.C) ...
             + prev_pos.nP*(cur.P - prev.P) ...
             + prev_pos.nF*(cur.F - prev.F);          % future = mark-to-market
        PnL_unh_cum = PnL_unh_cum + dV;
        PnL_hed_cum = PnL_hed_cum + (dV - dHed);

        % --- Rebalance: greeks at t_k -> new hedge -----------------------
        gC = greeks_vanilla_AB('call',   state_k.params_AB, state_k.scale_factor(2), state_k.mkt, bumps);
        gP = greeks_vanilla_AB('put',    state_k.params_AB, state_k.scale_factor(2), state_k.mkt, bumps);
        gF = greeks_vanilla_AB('future', state_k.params_AB, state_k.scale_factor(2), state_k.mkt, bumps);
        [pos_k, resG] = build_hedge_AB(cur.gEx, gC, gP, gF);

        trade = struct('nC', pos_k.nC - prev_pos.nC, ...
                       'nP', pos_k.nP - prev_pos.nP, ...
                       'nF', pos_k.nF - prev_pos.nF);
        cost_step = hedging_cost(trade, cur.F, costRule);
        cost_cum  = cost_cum + cost_step;

        rows{kk} = {snap, cur.Vex, PnL_unh_cum, PnL_hed_cum, cost_step, cost_cum, ...
                    resG.delta, resG.gamma, resG.vega};

        % roll forward
        prev     = cur;
        prev_pos = pos_k;
    end

    res = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'Date','V_exotic','PnL_unhedged','PnL_hedged', ...
         'cost_step','cost_cum','res_delta','res_gamma','res_vega'});
end

% =========================================================================
function m = mark_state(type, state, mc, bumps)
% Price the exotic (+ its greeks) and the two hedge vanillas + future level
% under a given calibrated `state`. Returns one struct of marks.
    gEx = greeks_exotic_AB(type, state.params_AB, state.scale_factor, state.mkt, mc, bumps);
    gC  = greeks_vanilla_AB('call', state.params_AB, state.scale_factor(2), state.mkt, []);
    gP  = greeks_vanilla_AB('put',  state.params_AB, state.scale_factor(2), state.mkt, []);
    m = struct('Vex', gEx.price, 'gEx', gEx, ...
               'C', gC.price, 'P', gP.price, 'F', state.mkt.forward);
end

% =========================================================================
function state = recalibrate_AB(snapDate, dataOpts, contract)
% Re-run EX1+EX2 on `snapDate` and repack the AB state for the SAME exotic.
%
%   Maturities are matched by EXPIRY DATE (contract.E1, contract.E2), not by
%   index: run_ex2 -> readData drops expired maturities and sorts by expiry,
%   so the positional index of a fixed calendar expiry can shift between
%   snapshots. The contract strikes (K1,K2,Kc,Kp) are FIXED terms and are
%   carried over unchanged; only forward / discount / scale are refreshed
%   (this is what makes theta show up in the hedged P&L).
    opts = struct('callpath',   dataOpts.callpath, ...
                  'putpath',     dataOpts.putpath, ...
                  'expiryFile',  dataOpts.expiryFile, ...
                  'valueDate',   snapDate, ...
                  'verbose',     false, ...
                  'plot',        false);
    [params, market] = run_ex2(opts);

    iT1 = match_maturity(market.expiries, contract.E1);
    iT2 = match_maturity(market.expiries, contract.E2);

    I0  = I0_AB(0, params.AB);
    sig = market.sigma_ATM / I0;
    scale_factor = [sig(iT1)*sqrt(market.yf(iT1)), sig(iT2)*sqrt(market.yf(iT2))];

    mkt = struct('forward', market.forward(iT2), ...
                 'K1', contract.K1, 'K2', contract.K2, ...
                 'Kc', contract.Kc, 'Kp', contract.Kp, ...
                 'df', [market.discount_factor(iT1), market.discount_factor(iT2)]);

    state = struct('params_AB', params.AB, 'scale_factor', scale_factor, 'mkt', mkt);
end

% =========================================================================
function idx = match_maturity(expiries, targetDate)
% Index of the row in `expiries` whose calendar date equals targetDate.
    [tol, idx] = min(abs(days(expiries - targetDate)));
    if tol > 0.5
        error('hedge_backtest:expiryNotFound', ...
            'Target expiry %s not present at this snapshot (closest off by %.1f days).', ...
            string(targetDate,'yyyy-MM-dd'), tol);
    end
end
