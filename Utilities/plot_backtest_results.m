function plot_backtest_results(LA_results_es6)
% PLOT_STRATEGY_COMPARISON Compares the 6 hedging strategies on common metrics.
%
% INPUTS:
%   LA_results_es6 - Struct containing results for test1...test6.

    % Extract available tests dynamically
    test_names = fieldnames(LA_results_es6);
    n_tests = length(test_names);
    
    % Initialize data vectors
    total_net      = zeros(1, n_tests);
    total_unhedged = zeros(1, n_tests);
    std_net        = zeros(1, n_tests);
    std_unhedged   = zeros(1, n_tests);
    cost_t0_arr    = zeros(1, n_tests);
    cost_rebal_arr = zeros(1, n_tests);
    test_labels    = cell(1, n_tests);
    
    % Custom descriptive labels (split onto two lines for clean plotting)
    custom_labels = {
        sprintf('Long Chooser 1'), ...
        sprintf('Long Chooser 2'), ...
        sprintf('Long CoC, PoP 1'), ...
        sprintf('Long CoC, PoP 2'), ...
        sprintf('Long CoC, Short Chooser 1'), ...
        sprintf('Long CoC, Short Chooser 2')
    };
    
    for i = 1:n_tests
        test_id = test_names{i};
        bt = LA_results_es6.(test_id).Backtest;
        
        if isempty(bt)
            continue;
        end
        
        total_net(i)      = bt.Total_Net_PnL;
        total_unhedged(i) = bt.Total_Unhedged_PnL;
        std_net(i)        = bt.Std_Hedged;
        std_unhedged(i)   = bt.Std_Unhedged;
        cost_t0_arr(i)    = bt.Cost_t0;
        cost_rebal_arr(i) = sum(bt.Costs);
        
        % Assign custom label if available, otherwise fallback to generic
        if i <= length(custom_labels)
            test_labels{i} = custom_labels{i};
        else
            test_labels{i} = sprintf('TEST %d', i);
        end
    end

    % =====================================================================
    % PLOT 1: FINAL TOTAL P&L COMPARISON
    % =====================================================================
    figure('Name', 'Strategy Comparison: Total P&L', 'Color', 'white', 'Position', [100, 100, 1000, 450]);
    
    bar_data = [total_unhedged', total_net'];
    b1 = bar(bar_data, 'grouped');
    
    b1(1).FaceColor = [0.8 0.2 0.2];
    b1(2).FaceColor = [0.2 0.6 0.2];
    
    set(gca, 'XTickLabel', test_labels, 'FontSize', 9);
    ylabel('Total P&L (EUR)', 'FontWeight', 'bold');
    title('Final Total P&L: Unhedged vs Hedged Across Strategies', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Unhedged P&L', 'Hedged Net P&L (After Costs)', 'Location', 'best');
    grid on;

    
    % =====================================================================
    % PLOT 2: TOTAL COSTS BREAKDOWN (STACKED)
    % =====================================================================
    figure('Name', 'Strategy Comparison: Costs Breakdown', 'Color', 'white', 'Position', [200, 200, 1000, 450]);
    
    cost_data = [cost_t0_arr', cost_rebal_arr'];
    b3 = bar(cost_data, 'stacked');
    
    b3(1).FaceColor = [0.3 0.3 0.3]; 
    b3(2).FaceColor = [0.9 0.6 0.1]; 
    
    set(gca, 'XTickLabel', test_labels, 'FontSize', 9);
    ylabel('Costs (EUR)', 'FontWeight', 'bold');
    title('Hedging Costs Breakdown: Initial Setup vs Dynamic Rebalancing', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Initial Setup Cost (t_0)', 'Dynamic Rebalancing Costs', 'Location', 'best');
    grid on;
    
    % Add total cost labels on top of the stacked bars
    for i = 1:n_tests
        total_cost = cost_t0_arr(i) + cost_rebal_arr(i);
        text(i, total_cost, sprintf(' %.0f €', total_cost), ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10, 'FontWeight', 'bold');
    end
end