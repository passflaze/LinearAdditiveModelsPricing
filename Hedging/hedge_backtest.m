function res = hedge_backtest(type, state0, positions0, tuesdays, dataOpts, mc, bumps, costRule)
% HEDGE_BACKTEST  Static-with-rebalancing hedge check over the next Tuesdays.
%
%   Builds the hedge at t0, then on each Tuesday re-observes the market,
%   RE-CALIBRATES the AB model, reprices the exotic and the hedge book,
%   rebalances the hedge and accumulates the bid-ask cost. Reports the P&L
%   of the unhedged vs hedged book and the residual greeks.
%
% INPUTS:
%   type       - (string) 'CoC' | 'PoP' | 'Chooser'
%   state0     - (struct) inception state at t0:
%                  .params_AB    [k; eta]
%                  .scale_factor [scale_t1, scale_t2]
%                  .mkt          (forward, K1, K2, df, Kc, Kp)   see price_exotic_AB
%   positions0 - (struct) .nC .nP .nF  hedge set at t0 (from build_hedge_AB)
%   tuesdays   - (Nx1 datetime) the value dates to test (next two Tuesdays)
%   dataOpts   - (struct) .callpath .putpath .expiryFile  (for readData)
%                          .iT1 .iT2  maturity indices of the (T1,T2) window
%   mc         - (struct) MC settings  (see price_exotic_AB)
%   bumps      - (struct) .dF .dSig
%   costRule   - (struct) .fut_bp .opt_bp
%
% OUTPUT:
%   res - (table) one row per date with columns:
%         Date, V_exotic, V_hedge, PnL_unhedged, PnL_hedged,
%         cost_step, cost_cum, res_delta, res_gamma, res_vega

    % --- t0 quantities ---------------------------------------------------
    V_exotic_0 = price_exotic_AB(type, state0.params_AB, state0.scale_factor, state0.mkt, mc);
    V_hedge_0  = hedge_value(state0, positions0, mc);
    cost_cum   = hedging_cost(positions0, state0.mkt.forward, costRule);   % inception cost

    prev_pos   = positions0;
    nT         = numel(tuesdays);
    rows       = cell(nT, 1);

    for kk = 1:nT
        snap = tuesdays(kk);

        % --- Re-observe market & re-calibrate AB at this Tuesday ----------
        state_k = recalibrate_AB(snap, dataOpts);   % see local function (TODO)

        % --- Reprice exotic and current hedge book -----------------------
        V_exotic_k = price_exotic_AB(type, state_k.params_AB, state_k.scale_factor, state_k.mkt, mc);
        V_hedge_k  = hedge_value(state_k, prev_pos, mc);    % OLD positions, NEW market

        % --- Greeks at Tuesday_k -> new hedge (rebalance) ----------------
        gEx = greeks_exotic_AB(type, state_k.params_AB, state_k.scale_factor, state_k.mkt, mc, bumps);
        gC  = greeks_vanilla_AB('call',   state_k.params_AB, state_k.scale_factor(2), state_k.mkt, bumps);
        gP  = greeks_vanilla_AB('put',    state_k.params_AB, state_k.scale_factor(2), state_k.mkt, bumps);
        gF  = greeks_vanilla_AB('future', state_k.params_AB, state_k.scale_factor(2), state_k.mkt, bumps);
        [pos_k, resG] = build_hedge_AB(gEx, gC, gP, gF);

        % --- Rebalancing trade & cost ------------------------------------
        trade = struct('nC', pos_k.nC - prev_pos.nC, ...
                       'nP', pos_k.nP - prev_pos.nP, ...
                       'nF', pos_k.nF - prev_pos.nF);
        cost_step = hedging_cost(trade, state_k.mkt.forward, costRule);
        cost_cum  = cost_cum + cost_step;

        % --- P&L of the book (long exotic, short hedge) ------------------
        %   unhedged : just the exotic mark-to-market move
        %   hedged   : exotic move minus hedge move (hedge held over the step)
        PnL_unhedged = V_exotic_k - V_exotic_0;
        PnL_hedged   = (V_exotic_k - V_exotic_0) - (V_hedge_k - V_hedge_0);

        rows{kk} = {snap, V_exotic_k, V_hedge_k, PnL_unhedged, PnL_hedged, ...
                    cost_step, cost_cum, resG.delta, resG.gamma, resG.vega};

        % roll forward: the rebalanced book becomes the reference for next step
        prev_pos   = pos_k;
        V_exotic_0 = V_exotic_k;
        V_hedge_0  = hedge_value(state_k, pos_k, mc);
    end

    res = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'Date','V_exotic','V_hedge','PnL_unhedged','PnL_hedged', ...
         'cost_step','cost_cum','res_delta','res_gamma','res_vega'});
end

% =========================================================================
function V = hedge_value(state, pos, mc) %#ok<INUSD>
% Mark-to-market value of the hedge book (call + put + future) under `state`.
    gC = greeks_vanilla_AB('call',   state.params_AB, state.scale_factor(2), state.mkt, []);
    gP = greeks_vanilla_AB('put',    state.params_AB, state.scale_factor(2), state.mkt, []);
    % Future MtM relative to its entry level is handled in P&L; here use spot F
    % as the per-unit value (delta-1 instrument). TODO(Persona B): if you track
    % the future entry price, value it as (F_now - F_entry); otherwise the
    % linear future term cancels in the P&L differences below.
    V = pos.nC * gC.price + pos.nP * gP.price + pos.nF * state.mkt.forward;
end

% =========================================================================
function state = recalibrate_AB(snapDate, dataOpts) %#ok<STOUT,INUSD>
% TODO(Persona B): re-run EX1+EX2 on `snapDate` and pack the AB state.
%   1) [strikes,calls,puts,expiries] = readData(dataOpts.callpath, ...
%          dataOpts.putpath, snapDate, dataOpts.expiryFile);
%   2) [params, market] = calibrate_surface(struct('valueDate',snapDate, ...));
%          (or run_ex2 with the same opts used in run_project2A, verbose=false)
%   3) Recompute the AB scale on the (iT1,iT2) window, exactly as run_ex4:
%          I0  = I0_AB(0, params.AB);
%          sig = market.sigma_ATM / I0;
%          scale_factor = [sig(iT1)*sqrt(market.yf(iT1)), sig(iT2)*sqrt(market.yf(iT2))];
%   4) Build mkt: forward = market.forward(iT2); K2 = forward; K1 unchanged;
%          df = [market.discount_factor(iT1), market.discount_factor(iT2)];
%          Kc, Kp = strikes of the hedge vanillas (keep them fixed from t0).
%   Return state = struct('params_AB',..,'scale_factor',..,'mkt',..).
    error('hedge_backtest:recalibrate_AB:notImplemented', ...
        'recalibrate_AB is a TODO — wire it to readData + calibrate_surface.');
end
