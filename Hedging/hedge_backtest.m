function res = hedge_backtest(state_t0, opts_array, prc_params, hedge_rules, greeks)
% HEDGE_BACKTEST  Static vs Dynamic backtest of a generic vanilla-basket hedge.
%
% INPUTS:
%   state_t0    - (struct) .V_exotic, .V_hedge, .w (Nx1 weights), .cost_t0
%   opts_array  - (struct array) options for run_ex2 (contains valueDate)
%   prc_params  - (struct) pricing setup, including .basket (Nx1 specs)
%   hedge_rules - (struct) .is_dynamic, .tol_delta/gamma/vega, .costRule
%   greeks      - (cell array) targeted Greeks, e.g. {'Delta','Gamma','Vega'}
%
% OUTPUT:
%   res         - (struct) PnL, costs and residual Greeks history.

    N_steps = length(opts_array);

    % Relative rebalancing tolerance: rebalance a Greek when the residual
    % exceeds rel_tol * |exotic Greek| (hedge error > rel_tol of gross exposure).
    if ~isfield(hedge_rules, 'rel_tol'), hedge_rules.rel_tol = 0.05; end

    % Initialize results arrays
    res.PnL         = zeros(N_steps, 1);   % hedged book (exotic + hedge)
    res.PnL_exotic  = zeros(N_steps, 1);   % UNHEDGED (exotic only)
    res.PnL_hedge   = zeros(N_steps, 1);   % hedge leg only
    res.Costs       = zeros(N_steps, 1);
    res.Net_PnL     = zeros(N_steps, 1);
    res.Book_Delta  = zeros(N_steps, 1);
    res.Book_Gamma  = zeros(N_steps, 1);
    res.Book_Vega   = zeros(N_steps, 1);

    % Day-0 hedging slippage (cost of opening the static hedge).
    if isfield(state_t0, 'cost_t0'), t0_cost = state_t0.cost_t0; else, t0_cost = 0; end

    % Track previous step's state
    old_V_exotic = state_t0.V_exotic;
    old_V_hedge  = state_t0.V_hedge;
    w_hedge      = state_t0.w(:);       % (N x 1) basket weights
    basket       = prc_params.basket;

    fprintf('\n--- INITIAL PORTFOLIO STATE (t0) ---\n');
    fprintf('  Exotic Portfolio Value : %+.4e\n', state_t0.V_exotic);
    fprintf('  Hedge Portfolio Value  : %+.4e\n', state_t0.V_hedge);

    for k = 1:N_steps
        current_date = opts_array(k).valueDate;

        fprintf('\n======================================================\n');
        fprintf(' EVALUATING TIME STEP %d: %s\n', k, datestr(current_date));
        fprintf('======================================================\n');

        % --- Market recalibration ---
        opts_current = opts_array(k);
        [params, market] = run_ex2(opts_current);
        iT1 = prc_params.iT1;
        iT2 = prc_params.iT2;

        I0_AB_val  = I0_AB(0, params.AB);
        sigma_t_AB = market.sigma_ATM / I0_AB_val;

        scale_factor_exotic  = [sigma_t_AB(iT1)*sqrt(market.yf(iT1)), sigma_t_AB(iT2)*sqrt(market.yf(iT2))];
        scale_factor_vanilla = sigma_t_AB .* sqrt(market.yf);

        % Forward F(t0,T2) used by the exotic pricer / vega recalibration.
        prc_params.params_hedge.forward = market.forward(iT2);

        % --- Reprice exotic portfolio ---
        nE = numel(prc_params.exotics);
        G_exotic = repmat(struct('price',0,'delta',0,'gamma',0,'vega',0), nE, 1);

        for ee = 1:nE
            G_exotic(ee) = greeks_exotic_AB(prc_params.exotics{ee}, params.AB, scale_factor_exotic, market, prc_params.params_hedge, prc_params.mc, prc_params.bumps);
        end

        w_exo = prc_params.w_exotics;
        port_exotic = struct('price',0, 'delta',0, 'gamma',0, 'vega',0);
        for ee = 1:nE
            port_exotic.price = port_exotic.price + w_exo(ee) * G_exotic(ee).price;
            port_exotic.delta = port_exotic.delta + w_exo(ee) * G_exotic(ee).delta;
            port_exotic.gamma = port_exotic.gamma + w_exo(ee) * G_exotic(ee).gamma;
            port_exotic.vega  = port_exotic.vega  + w_exo(ee) * G_exotic(ee).vega;
        end
        new_V_exotic = port_exotic.price;

        % --- Reprice hedging basket ---
        instr   = price_hedge_basket(basket, params.AB, scale_factor_vanilla, market, ...
                                     prc_params.params_hedge, prc_params.mc, prc_params.bumps);
        prices  = [instr.price]';
        i_delta = [instr.delta]';
        i_gamma = [instr.gamma]';
        i_vega  = [instr.vega]';

        % Mark-to-market value of the existing hedge basket.
        current_V_old_hedge = w_hedge' * prices;

        % --- Gross P&L of the hedged book ---
        fprintf('  Exotic Portfolio MTM : %+.4e\n', new_V_exotic);
        fprintf('  Hedge Portfolio MTM  : %+.4e\n', current_V_old_hedge);

        dPnL_exotic = new_V_exotic - old_V_exotic;
        dPnL_hedge  = current_V_old_hedge - old_V_hedge;

        res.PnL_exotic(k) = dPnL_exotic;          % unhedged
        res.PnL_hedge(k)  = dPnL_hedge;
        res.PnL(k)        = dPnL_exotic + dPnL_hedge;
        fprintf('  Unhedged P&L (exotic): %+.4e\n', dPnL_exotic);
        fprintf('  Gross P&L over period: %+.4e\n', res.PnL(k));

        % --- Dynamic hedging check (rebalancing) ---
        book = struct( ...
            'delta', port_exotic.delta + w_hedge' * i_delta, ...
            'gamma', port_exotic.gamma + w_hedge' * i_gamma, ...
            'vega',  port_exotic.vega  + w_hedge' * i_vega);

        res.Book_Delta(k) = book.delta;
        res.Book_Gamma(k) = book.gamma;
        res.Book_Vega(k)  = book.vega;
        fprintf('  Residual Book Delta  : %+.4e\n', book.delta);
        fprintf('  Residual Book Gamma  : %+.4e\n', book.gamma);
        fprintf('  Residual Book Vega   : %+.4e\n', book.vega);

        transaction_cost = 0;

        % Trigger if ANY targeted Greek's residual exceeds rel_tol of the
        % gross exotic exposure to that Greek.
        breach = false;
        for gi = 1:numel(greeks)
            gn  = lower(greeks{gi});
            tol = hedge_rules.rel_tol * abs(port_exotic.(gn));
            if abs(book.(gn)) > tol
                breach = true;
                fprintf('  %s residual %+.2e exceeds %.0f%% of exposure (tol %.2e).\n', ...
                        greeks{gi}, book.(gn), 100*hedge_rules.rel_tol, tol);
            end
        end

        if hedge_rules.is_dynamic && breach
            fprintf('  Rebalancing...\n');

            [w_new, ~] = build_hedge(port_exotic, instr, greeks);

            trade_qty = w_new - w_hedge;

            costs = hedging_cost(trade_qty, instr, hedge_rules.costRule);
            transaction_cost = costs.slippage;

            fprintf('  Rebalanced! Slippage: $%+.4e | Premium: $%+.4e | Total Cash Flow: $%+.4e\n', ...
                    costs.slippage, costs.premium, costs.total);

            w_hedge = w_new;
        end

        % --- Store results and advance state ---
        res.Costs(k)   = transaction_cost;
        res.Net_PnL(k) = res.PnL(k) - transaction_cost;

        old_V_exotic = new_V_exotic;
        old_V_hedge  = w_hedge' * prices;   % post-rebalance weights at current prices
    end

    % Fold the day-0 hedge opening cost into the totals.
    res.Cost_t0       = t0_cost;
    total_costs       = sum(res.Costs) + t0_cost;
    total_net_pnl     = sum(res.Net_PnL) - t0_cost;
    res.Total_Net_PnL = total_net_pnl;

    % Unhedged benchmark: just holding the exotic portfolio (no hedge, no cost).
    res.Total_Unhedged_PnL = sum(res.PnL_exotic);

    fprintf('\n=== BACKTEST COMPLETE ===\n');
    fprintf('Total Gross P&L     : %+.4e\n', sum(res.PnL));
    fprintf('Day-0 hedging cost  : %+.4e\n', t0_cost);
    fprintf('Rebalancing costs   : %+.4e\n', sum(res.Costs));
    fprintf('Total Costs         : %+.4e\n', total_costs);
    fprintf('Total Net P&L       : %+.4e\n', total_net_pnl);

    fprintf('\n--- HEDGED vs UNHEDGED ---\n');
    fprintf('Unhedged P&L (exotic only)   : %+.4e\n', res.Total_Unhedged_PnL);
    fprintf('Hedged Net P&L (after costs) : %+.4e\n', total_net_pnl);
    

end
