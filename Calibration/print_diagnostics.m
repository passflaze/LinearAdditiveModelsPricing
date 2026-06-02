function print_diagnostics(k_AB, eta_AB, fval_AB, alpha_MA, beta_MA, fval_MA, ...
        alpha_GL, beta_GL, fval_GL, M, dz, discount_factor, yf, sigma_ATM, ...
        moneyness_modified, c_mkt_calibration, expiries, nT, plot_flag)
% Post-calibration diagnostics and reporting (verbose path).
%
% Added plot_flag (optional, default false): if true, plots a comparison 
% of model vs market prices and errors in bps for the sampled options.

if nargin < 19
    plot_flag = false;
end

fprintf('STEP 5: Post-Calibration Diagnostics and Reporting...\n');
% Compute model prices using calibrated parameters
c_mod_AB = price_AB([k_AB, eta_AB], discount_factor, yf, sigma_ATM, moneyness_modified);
c_mod_MA = price_MA([alpha_MA, beta_MA], discount_factor, yf, sigma_ATM, moneyness_modified);
c_mod_GL = price_GL(alpha_GL, beta_GL, M, dz, discount_factor, sigma_ATM, yf, moneyness_modified);

% Clamp tiny negative artefacts from FFT/interp extrapolation to 0.
c_mod_AB(c_mod_AB < 0) = 0;
c_mod_MA(c_mod_MA < 0) = 0;
c_mod_GL(c_mod_GL < 0) = 0;

% Calculate residuals
res_AB = c_mkt_calibration - c_mod_AB;
res_MA = c_mkt_calibration - c_mod_MA;
res_GL = c_mkt_calibration - c_mod_GL;

% Print Parameter Summary Table
fprintf('\n=========================================================================\n');
fprintf('                     CALIBRATED PARAMETERS SUMMARY                       \n');
fprintf('=========================================================================\n');
fprintf('%-10s | %-12s | %-12s | %-15s\n', 'Model', 'Param 1', 'Param 2', 'Final SSE');
fprintf('-------------------------------------------------------------------------\n');
fprintf('%-10s | %-12.6f | %-12.6f | %-15.6g\n', 'AB (k,eta)', k_AB, eta_AB, fval_AB);
fprintf('%-10s | %-12.6f | %-12.6f | %-15.6g\n', 'MA (a,b)', alpha_MA, beta_MA, fval_MA);
fprintf('%-10s | %-12.6f | %-12.6f | %-15.6g\n', 'GL (a,b)', alpha_GL, beta_GL, fval_GL);

% Print RMSE per Maturity Table
fprintf('\n=========================================================================\n');
fprintf('                     ROOT MEAN SQUARE ERROR (RMSE)                       \n');
fprintf('=========================================================================\n');
fprintf('%-15s | %-12s | %-12s | %-12s\n', 'Expiry Date', 'RMSE AB', 'RMSE MA', 'RMSE GL');
fprintf('-------------------------------------------------------------------------\n');
for i = 1:nT
    r_AB_i = res_AB(i,:); r_AB_i = r_AB_i(~isnan(r_AB_i));
    r_MA_i = res_MA(i,:); r_MA_i = r_MA_i(~isnan(r_MA_i));
    r_GL_i = res_GL(i,:); r_GL_i = r_GL_i(~isnan(r_GL_i));
    fprintf('%-15s | %-12.6f | %-12.6f | %-12.6f\n', string(expiries(i),'yyyy-MM-dd'), ...
        sqrt(mean(r_AB_i.^2)), sqrt(mean(r_MA_i.^2)), sqrt(mean(r_GL_i.^2)));
end

% Print Sample Pricing Comparison Table (Fixed Maturity: Row 4)
fprintf('\n=============================================================================================================\n');
fprintf('                 SAMPLE PRICING COMPARISON (MARKET VS MODELS) - MATURITY 1 YEAR                              \n');
fprintf('=============================================================================================================\n');
fprintf('%-10s | %-10s | %-12s | %-12s | %-12s | %-12s | %-13s | %-13s | %-13s\n', 'Maturity', 'Moneyness', 'Market Px', 'AB Price', 'MA Price', 'GL Price', 'Diff AB (bps)', 'Diff MA (bps)', 'Diff GL (bps)');
fprintf('-------------------------------------------------------------------------------------------------------------\n');

% Fix row 4 for 1-year maturity
target_r = 4;
% Find valid columns (non-NaN) for the specified maturity
valid_cols = find(~isnan(c_mkt_calibration(target_r, :)));
num_samples = min(10, length(valid_cols));

% Select roughly evenly spaced moneyness samples from the valid columns
sample_col_indices = round(linspace(1, length(valid_cols), num_samples));
selected_cols = valid_cols(sample_col_indices);

% Initialize arrays for plotting if flag is true
if plot_flag
    x_labels   = strings(num_samples, 1);
    px_mkt     = zeros(num_samples, 1);
    px_mod_AB  = zeros(num_samples, 1);
    px_mod_MA  = zeros(num_samples, 1);
    px_mod_GL  = zeros(num_samples, 1);
    err_AB_bps = zeros(num_samples, 1);
    err_MA_bps = zeros(num_samples, 1);
    err_GL_bps = zeros(num_samples, 1);
end

for k = 1:num_samples
    c = selected_cols(k);
    
    % Calculate differences in basis points (Model - Market)
    diff_AB_bps = (c_mod_AB(target_r, c) - c_mkt_calibration(target_r, c)) * 10000;
    diff_MA_bps = (c_mod_MA(target_r, c) - c_mkt_calibration(target_r, c)) * 10000;
    diff_GL_bps = (c_mod_GL(target_r, c) - c_mkt_calibration(target_r, c)) * 10000;
    
    fprintf('%-10s | %-10.4f | %-12.6f | %-12.6f | %-12.6f | %-12.6f | %-13.2f | %-13.2f | %-13.2f\n', ...
        string(expiries(target_r), 'yyyy-MM'), ...
        moneyness_modified(target_r, c), ...
        c_mkt_calibration(target_r, c), ...
        c_mod_AB(target_r, c), ...
        c_mod_MA(target_r, c), ...
        c_mod_GL(target_r, c), ...
        diff_AB_bps, ...
        diff_MA_bps, ...
        diff_GL_bps);
        
    % Store data for plotting
    if plot_flag
        % Simplified label since maturity is fixed
        x_labels(k)   = sprintf('x=%.1f', moneyness_modified(target_r, c)); 
        px_mkt(k)     = c_mkt_calibration(target_r, c);
        px_mod_AB(k)  = c_mod_AB(target_r, c);
        px_mod_MA(k)  = c_mod_MA(target_r, c);
        px_mod_GL(k)  = c_mod_GL(target_r, c);
        err_AB_bps(k) = diff_AB_bps;
        err_MA_bps(k) = diff_MA_bps;
        err_GL_bps(k) = diff_GL_bps;
    end
end
fprintf('=============================================================================================================\n');

% =========================================================================
% OPTIONAL PLOT: MARKET VS MODELS & BPS ERRORS (FIXED MATURITY)
% =========================================================================
if plot_flag
    figure('Name', 'Diagnostics: Sample Pricing & Errors', 'Position', [100, 100, 1000, 800], 'Color', 'w');
    
    expiry_str = string(expiries(target_r), 'yyyy-MM-dd');
    
    % --- Subplot 1: Prices ---
    subplot(2, 1, 1);
    hold on; grid on; box on;
    plot(1:num_samples, px_mkt, 'ko-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', 'Market Price');
    plot(1:num_samples, px_mod_AB, 'r^-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'AB Price');
    plot(1:num_samples, px_mod_MA, 'bs-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'MA Price');
    plot(1:num_samples, px_mod_GL, 'g*-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'GL Price');
    
    title(sprintf('Pricing Comparison vs Market (Maturity: %s)', expiry_str), 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Option Price ($)', 'FontSize', 12);
    xticks(1:num_samples);
    xticklabels(x_labels);
    legend('Location', 'best', 'FontSize', 11);
    hold off;
    
    % --- Subplot 2: Errors in BPS ---
    subplot(2, 1, 2);
    hold on; grid on; box on;
    plot(1:num_samples, err_AB_bps, 'r^-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Error AB');
    plot(1:num_samples, err_MA_bps, 'bs-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Error MA');
    plot(1:num_samples, err_GL_bps, 'g*-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Error GL');
    
    % Add a zero-reference line
    yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
    
    title(sprintf('Pricing Error (Model - Market) in Bps (Maturity: %s)', expiry_str), 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Difference (bps)', 'FontSize', 12);
    xticks(1:num_samples);
    xticklabels(x_labels);
    legend('Location', 'best', 'FontSize', 11);
    hold off;
end

end