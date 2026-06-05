function plot_backtest_results(LA_results_es6)
% PLOT_BACKTEST_RESULTS  Publication-quality hedging backtest summary.
%
% Produces:
%   Figure 1 – Per-strategy P&L panels       (small multiples)
%   Figure 2 – Per-strategy cost panels      (small multiples)
%   Figure 3 – Rendered summary table        (manually drawn)
%   Console  – Formatted strategy comparison table

    % ── Design tokens ────────────────────────────────────────────────────
    % Colour
    C_uh  = [0.86 0.45 0.19];    % terracotta    – unhedged
    C_hd  = [0.20 0.40 0.68];    % oxford blue   – hedged
    C_c0  = [0.36 0.36 0.44];    % warm charcoal – setup cost
    C_cr  = [0.88 0.68 0.20];    % gold          – rebalancing
    DARK  = [0.10 0.10 0.14];    % near-black
    MID   = [0.50 0.50 0.54];    % medium gray
    LIGHT = [0.82 0.82 0.85];    % light gray
    PANEL = [0.973 0.973 0.978]; % panel background
    WHITE = [1.00 1.00 1.00];
    HDR   = [0.12 0.18 0.30];    % dark-navy header
    C_pos = [0.07 0.40 0.17];    % dark green (positive)
    C_neg = [0.58 0.10 0.10];    % dark red   (negative)

    % ── Extract test data ─────────────────────────────────────────────────
    test_names = fieldnames(LA_results_es6);
    n_tests    = numel(test_names);

    total_net      = zeros(1, n_tests);
    total_unhedged = zeros(1, n_tests);
    std_net        = zeros(1, n_tests);
    std_unhedged   = zeros(1, n_tests);
    cost_t0_arr    = zeros(1, n_tests);
    cost_rebal_arr = zeros(1, n_tests);
    test_labels    = cell(1, n_tests);

    custom_labels = { ...
        'Long Chooser 1',            ...
        'Long Chooser 2',            ...
        'Long CoC, PoP 1',           ...
        'Long CoC, PoP 2',           ...
        'Long CoC, Short Chooser 1', ...
        'Long CoC, Short Chooser 2'};

    for i = 1:n_tests
        bt = LA_results_es6.(test_names{i}).Backtest;
        if isempty(bt); continue; end
        total_net(i)      = bt.Total_Net_PnL;
        total_unhedged(i) = bt.Total_Unhedged_PnL;
        std_net(i)        = bt.Std_Hedged;
        std_unhedged(i)   = bt.Std_Unhedged;
        cost_t0_arr(i)    = bt.Cost_t0;
        cost_rebal_arr(i) = sum(bt.Costs);
        if i <= numel(custom_labels)
            test_labels{i} = custom_labels{i};
        else
            test_labels{i} = sprintf('Test %d', i);
        end
    end

    total_cost = cost_t0_arr + cost_rebal_arr;
    var_red    = (std_unhedged - std_net) ./ max(std_unhedged, eps) * 100;
    delta_pnl  = total_net - total_unhedged;

    n_gc = min(n_tests, 3);           % grid cols
    n_gr = ceil(n_tests / n_gc);      % grid rows
    w_s  = 0.30;                      % bar centre offset from x=1
    bw   = 0.38;                      % bar width

    % ══════════════════════════════════════════════════════════════════════
    % FIGURE 1 – Per-Strategy P&L
    % ══════════════════════════════════════════════════════════════════════
    fig1 = figure('Name', 'Per-Strategy P&L', 'Color', WHITE, ...
                  'Position', [60 100 n_gc*340 n_gr*325]);
    tl1  = tiledlayout(fig1, n_gr, n_gc, 'TileSpacing', 'tight', 'Padding', 'compact');
    sgtitle(tl1, 'P&L per Strategy — Unhedged vs. Hedged', ...
            'FontSize', 14, 'FontWeight', 'bold', 'Color', DARK);

    for i = 1:n_tests
        ax = nexttile(tl1);
        set(ax, 'Color', PANEL);
        hold(ax, 'on');

        % Bars with subtle dark edge
        bar(ax, 1-w_s, total_unhedged(i), bw, ...
            'FaceColor', C_uh, 'EdgeColor', C_uh*0.72, 'LineWidth', 0.5);
        bar(ax, 1+w_s, total_net(i),      bw, ...
            'FaceColor', C_hd, 'EdgeColor', C_hd*0.72, 'LineWidth', 0.5);

        % Error bars – no marker dot, mid-gray
        errorbar(ax, 1-w_s, total_unhedged(i), std_unhedged(i), ...
                 'Color', MID, 'LineWidth', 0.85, 'CapSize', 3, 'Marker', 'none');
        errorbar(ax, 1+w_s, total_net(i),      std_net(i), ...
                 'Color', MID, 'LineWidth', 0.85, 'CapSize', 3, 'Marker', 'none');

        % Zero reference line
        yline(ax, 0, 'Color', LIGHT, 'LineWidth', 1.1);

        % Value labels, placed beyond error caps, colour-coded ±
        max_y_i = max(abs([total_net(i), total_unhedged(i), std_net(i), std_unhedged(i)])) + eps;
        nudge   = 0.06 * max_y_i;

        y_u  = total_unhedged(i) + sign(total_unhedged(i)+eps)*(std_unhedged(i)+nudge);
        va_u = 'bottom'; if total_unhedged(i) < 0, va_u = 'top'; end
        lc_u = C_uh*0.68; if total_unhedged(i) < 0, lc_u = C_neg; end
        text(ax, 1-w_s, y_u, sprintf('%+.0f', total_unhedged(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', va_u, ...
             'FontSize', 7.5, 'FontWeight', 'bold', 'Color', lc_u);

        y_h  = total_net(i) + sign(total_net(i)+eps)*(std_net(i)+nudge);
        va_h = 'bottom'; if total_net(i) < 0, va_h = 'top'; end
        lc_h = C_hd*0.72; if total_net(i) < 0, lc_h = C_neg; end
        text(ax, 1+w_s, y_h, sprintf('%+.0f', total_net(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', va_h, ...
             'FontSize', 7.5, 'FontWeight', 'bold', 'Color', lc_h);

        % Variance-reduction badge at top centre
        red_i   = (std_unhedged(i) - std_net(i)) / max(std_unhedged(i), eps) * 100;
        ax_ylim = ylim(ax);
        text(ax, 1.0, ax_ylim(2) - 0.03*(ax_ylim(2)-ax_ylim(1)), ...
             sprintf('\\downarrow %.0f%%', red_i), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
             'FontSize', 8, 'Color', MID, 'Interpreter', 'tex');

        % Axis styling
        set(ax, ...
            'XTick',      [1-w_s, 1+w_s], ...
            'XTickLabel', {'Unhedged', 'Hedged'}, ...
            'FontSize',   8.5, ...
            'Color',      PANEL, 'GridColor', LIGHT, 'GridAlpha', 1, ...
            'Box',        'off', 'TickDir', 'out', ...
            'XGrid',      'off', 'YGrid', 'on', ...
            'TickLength', [0.025 0.025], ...
            'LineWidth',  0.55, 'XColor', MID, 'YColor', MID);
        xlim(ax, [1 - w_s - bw/2 - 0.10,  1 + w_s + bw/2 + 0.10]);
        ax.YAxis.Exponent = 0;
        title(ax, test_labels{i}, 'FontSize', 10.5, 'FontWeight', 'bold', 'Color', DARK);
        ylabel(ax, 'P&L (EUR)', 'FontSize', 8, 'Color', MID);
    end

    % ══════════════════════════════════════════════════════════════════════
    % FIGURE 2 – Per-Strategy Cost Panels
    % ══════════════════════════════════════════════════════════════════════
    fig2 = figure('Name', 'Per-Strategy Costs', 'Color', WHITE, ...
                  'Position', [200 130 n_gc*340 n_gr*310]);
    tl2  = tiledlayout(fig2, n_gr, n_gc, 'TileSpacing', 'tight', 'Padding', 'compact');
    sgtitle(tl2, 'Hedging Costs per Strategy', ...
            'FontSize', 14, 'FontWeight', 'bold', 'Color', DARK);

    for i = 1:n_tests
        ax = nexttile(tl2);
        set(ax, 'Color', PANEL);
        hold(ax, 'on');

        bc_i = bar(ax, 1, [cost_t0_arr(i), cost_rebal_arr(i)], 0.48, ...
                   'stacked', 'EdgeColor', 'none');
        bc_i(1).FaceColor = C_c0;
        bc_i(2).FaceColor = C_cr;

        total_i = cost_t0_arr(i) + cost_rebal_arr(i);
        nudge_i = 0.06 * (total_i + eps);

        % Total label above bar
        text(ax, 1, total_i + nudge_i, sprintf('%.0f EUR', total_i), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 9.5, 'FontWeight', 'bold', 'Color', DARK);

        % Segment labels (only if segment big enough to hold them)
        thresh = 0.09 * (total_i + eps);
        if cost_t0_arr(i) > thresh
            text(ax, 1, cost_t0_arr(i)/2, sprintf('%.0f', cost_t0_arr(i)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 8, 'Color', WHITE, 'FontWeight', 'bold');
        end
        if cost_rebal_arr(i) > thresh
            text(ax, 1, cost_t0_arr(i) + cost_rebal_arr(i)/2, ...
                 sprintf('%.0f', cost_rebal_arr(i)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 8, 'Color', WHITE, 'FontWeight', 'bold');
        end

        set(ax, 'XTick', [], 'FontSize', 8.5, ...
            'Color', PANEL, 'GridColor', LIGHT, 'GridAlpha', 1, ...
            'Box', 'off', 'TickDir', 'out', 'XGrid', 'off', 'YGrid', 'on', ...
            'TickLength', [0.025 0.025], ...
            'LineWidth', 0.55, 'XColor', MID, 'YColor', MID);
        xlim(ax, [0.5, 1.5]);
        ax.YAxis.Exponent = 0;
        title(ax, test_labels{i}, 'FontSize', 10.5, 'FontWeight', 'bold', 'Color', DARK);
        ylabel(ax, 'Cost (EUR)', 'FontSize', 8, 'Color', MID);

        if i == 1
            lg = legend(ax, bc_i, {'Setup (t_0)', 'Rebalancing'}, ...
                        'Location', 'north', 'FontSize', 8.5, 'Box', 'off');
            lg.TextColor = DARK;
        end
    end

    % ══════════════════════════════════════════════════════════════════════
    % CONSOLE TABLE
    % ══════════════════════════════════════════════════════════════════════
    sep_t = repmat('-', 1, 110);
    sep_b = repmat('=', 1, 110);
    fprintf('\n%s\n', sep_b);
    fprintf('  HEDGING BACKTEST — STRATEGY COMPARISON\n');
    fprintf('%s\n', sep_b);
    fprintf('  %-28s  %10s  %10s  %10s  %10s  %10s  %8s  %10s\n', ...
            'Strategy', 'Uhgd P&L', 'Hgd P&L', 'Delta P&L', ...
            'sig Unhdg', 'sig Hedgd', 'VarRed%', 'TotCost');
    fprintf('%s\n', sep_t);
    for i = 1:n_tests
        fprintf('  %-28s  %+10.1f  %+10.1f  %+10.1f  %10.1f  %10.1f  %7.1f%%  %10.1f\n', ...
                test_labels{i}, total_unhedged(i), total_net(i), delta_pnl(i), ...
                std_unhedged(i), std_net(i), var_red(i), total_cost(i));
    end
    fprintf('%s\n\n', sep_b);

    % ══════════════════════════════════════════════════════════════════════
    % FIGURE 3 – Manually rendered summary table
    % ══════════════════════════════════════════════════════════════════════
    px_t = 56;   % title height
    px_h = 40;   % header row
    px_r = 38;   % data row
    px_b = 18;   % bottom padding
    fw   = 1150; % figure width (px)
    fh   = px_t + px_h + n_tests * px_r + px_b;

    fig3 = figure('Name', 'Backtest Summary Table', 'Color', WHITE, ...
                  'Position', [350 200 fw fh]);
    ax_t = axes(fig3, 'Position', [0 0 1 1], ...
                'XLim', [0 1], 'YLim', [0 1], 'Visible', 'off');
    hold(ax_t, 'on');

    % ── Helper: pixel-from-top → normalised y ────────────────────────────
    n2y = @(p) 1 - p/fh;

    % ── Table title ───────────────────────────────────────────────────────
    text(ax_t, 0.5, n2y(px_t/2), 'Hedging Backtest — Strategy Comparison', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 13, 'FontWeight', 'bold', 'Color', DARK);

    % ── Column layout ─────────────────────────────────────────────────────
    col_x     = [0.000, 0.170, 0.260, 0.350, 0.440, 0.530, 0.618, 0.700, 0.778, 0.876, 1.000];
    col_names = {'Strategy',   'Unhdg P&L', 'Hgd P&L', 'Delta P&L', ...
                 'sig Unhdg',  'sig Hedgd', 'VarRed %', ...
                 'Cost t0',    'Rebal',     'Total Cost'};

    % ── Header row ────────────────────────────────────────────────────────
    hdr_y1 = n2y(px_t);
    hdr_y0 = n2y(px_t + px_h);
    patch(ax_t, [0 1 1 0], [hdr_y0 hdr_y0 hdr_y1 hdr_y1], HDR, 'EdgeColor', 'none');

    for j = 1:10
        xc = (col_x(j) + col_x(j+1)) / 2;
        yc = (hdr_y0 + hdr_y1) / 2;
        text(ax_t, xc, yc, col_names{j}, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', 8.5, 'FontWeight', 'bold', 'Color', WHITE);
    end

    % ── Data rows ─────────────────────────────────────────────────────────
    STRIPE = [0.955 0.958 0.966];
    rh_n   = px_r / fh;   % row height in normalised units

    for i = 1:n_tests
        r_y1 = hdr_y0 - (i-1)*rh_n;
        r_y0 = r_y1 - rh_n;
        r_yc = (r_y0 + r_y1) / 2;

        row_bg = WHITE; if mod(i,2) == 0, row_bg = STRIPE; end
        patch(ax_t, [0 1 1 0], [r_y0 r_y0 r_y1 r_y1], row_bg, 'EdgeColor', LIGHT);

        vals = { test_labels{i}, ...
                 sprintf('%+.1f', total_unhedged(i)), ...
                 sprintf('%+.1f', total_net(i)),      ...
                 sprintf('%+.1f', delta_pnl(i)),      ...
                 sprintf('%.1f',  std_unhedged(i)),   ...
                 sprintf('%.1f',  std_net(i)),        ...
                 sprintf('%.1f%%', var_red(i)),       ...
                 sprintf('%.1f',  cost_t0_arr(i)),    ...
                 sprintf('%.1f',  cost_rebal_arr(i)), ...
                 sprintf('%.1f',  total_cost(i)) };

        % Value-based colour coding (cols 2-4 = P&L series, 7 = VarRed)
        pnl_sign  = [total_unhedged(i), total_net(i), delta_pnl(i)];
        vred_sign = var_red(i);

        for j = 1:10
            if j == 1
                xp = col_x(j) + 0.006; ha = 'left';
            else
                xp = col_x(j+1) - 0.006; ha = 'right';
            end

            tc = DARK;
            if j >= 2 && j <= 4
                v = pnl_sign(j-1);
                if v > 0, tc = C_pos; elseif v < 0, tc = C_neg; end
            elseif j == 7
                if vred_sign > 0, tc = C_pos; end
            end

            fw_str = 'normal';
            if j == 10, fw_str = 'bold'; end

            text(ax_t, xp, r_yc, vals{j}, ...
                 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle', ...
                 'FontSize', 8.5, 'Color', tc, 'FontWeight', fw_str);
        end
    end

    % ── Table borders ─────────────────────────────────────────────────────
    bot_y = hdr_y0 - n_tests * rh_n;

    % Bottom border
    line(ax_t, [0.005, 0.995], [bot_y, bot_y], 'Color', MID, 'LineWidth', 0.9);

    % Vertical column separators
    for j = 2:10
        line(ax_t, [col_x(j), col_x(j)], [bot_y, hdr_y1], 'Color', LIGHT, 'LineWidth', 0.5);
    end

    % Fix axes limits after all drawing
    xlim(ax_t, [0 1]);
    ylim(ax_t, [0 1]);

end
