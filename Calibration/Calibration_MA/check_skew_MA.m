function skew_report = check_skew_MA(alpha_MA, beta_MA, discount_factor, yf, ...
                                     sigma_ATM, moneyness_modified, ...
                                     c_mkt_calibration, expiries)
% CHECK_SKEW_MA  Post-calibration skew diagnostic for the Minimal Additive model.
%
%   Checks the MA fit via the Bachelier IV skew. The MA driver is an asymmetric
%   Laplace whose skewness depends only on beta/alpha (the single identifiable
%   parameter), and by separability the normalized smile w(chi)=sigma_imp/sigma_ATM
%   is maturity-independent in chi = (K-F)/(sigma_ATM*sqrt(t)). For each maturity
%   we invert prices to sigma_imp and fit w ~ a0 + a1*chi + a2*chi^2 on
%   |chi| <= CHI_WIN: ATM skew = a1, level = a0 (~1), curvature = 2*a2. The model
%   skew is measured identically on a dense chi grid.
%
%   INPUTS
%     alpha_MA, beta_MA   : calibrated MA parameters (scalars)
%     discount_factor, yf, sigma_ATM : (M x 1) per maturity
%     moneyness_modified  : (M x N) chi (NaN where no data)
%     c_mkt_calibration   : (M x N) market call prices (NaN where no data)
%     expiries            : (M x 1) datetime expiry dates (labels)
%   OUTPUT
%     skew_report : struct (chi_grid, w_model, skew/level/curv model & market,
%                   gamma1_dist analytic skewness, invariance_dev).
%   Reference: Baviera & Massaria (2026), Eq. (15) and separability.

    % --------------------------- settings --------------------------------
    CHI_WIN = 2.0;     % half-width of the ATM window (in chi) for the local fit
    NGRID   = 401;     % points on the dense model chi grid

    alpha_MA = alpha_MA(:); beta_MA = beta_MA(:);
    sigma_ATM = sigma_ATM(:); yf = yf(:); discount_factor = discount_factor(:);
    M = numel(yf);

    % ===================================================================== %
    %  1) ANALYTIC ASYMMETRIC-LAPLACE SKEWNESS (theory, distribution-level)  %
    % ===================================================================== %
    % gamma1 = 2*(1/b^3 - 1/a^3)/(1/a^2 + 1/b^2)^(3/2); scale-invariant in b/a.
    a = alpha_MA; b = beta_MA;
    gamma1_dist = 2*(1/b^3 - 1/a^3) / (1/a^2 + 1/b^2)^(3/2);

    % ===================================================================== %
    %  2) MARKET SKEW PER MATURITY (invert prices -> normalized IV -> fit)   %
    % ===================================================================== %
    skew_mkt  = nan(M, 1);
    level_mkt = nan(M, 1);
    curv_mkt  = nan(M, 1);

    % Invert the whole surface to Bachelier IV at once, then normalize per
    % maturity by sigma_ATM. Dollar moneyness F-K = -chi*sigma_ATM*sqrt(t).
    N           = size(moneyness_modified, 2);
    chi_mkt_all = moneyness_modified;
    sT          = sigma_ATM .* sqrt(yf);                   % (M x 1)
    fmk_all     = -moneyness_modified .* sT;               % (M x N)
    B_all       = repmat(discount_factor, 1, N);          % (M x N)
    t_all       = repmat(yf, 1, N);                        % (M x N)
    valid_all   = ~isnan(moneyness_modified) & ~isnan(c_mkt_calibration) & (c_mkt_calibration > 0);

    sig_imp_all          = nan(size(moneyness_modified));
    sig_imp_all(valid_all) = bachelier_iv(c_mkt_calibration(valid_all), B_all(valid_all), ...
                                          fmk_all(valid_all), t_all(valid_all));
    w_mkt_all = sig_imp_all ./ sigma_ATM;                  % (M x N) normalized smile

    % Local quadratic ATM fit per maturity (reads the precomputed smile).
    for i = 1:M
        v = ~isnan(chi_mkt_all(i, :)) & ~isnan(w_mkt_all(i, :));
        if nnz(v) < 3, continue; end
        [a0, a1, a2] = local_quad_fit(chi_mkt_all(i, v), w_mkt_all(i, v), CHI_WIN);
        level_mkt(i) = a0;
        skew_mkt(i)  = a1;
        curv_mkt(i)  = 2*a2;
    end

    % ===================================================================== %
    %  3) MODEL SKEW (MA smile on a dense chi grid, maturity-independent)    %
    % ===================================================================== %
    chi_grid = linspace(-CHI_WIN*1.5, CHI_WIN*1.5, NGRID);

    % Representative (median) maturity; smile is maturity-independent by separability.
    iref = max(1, round(M/2));

    w_model = ma_normalized_smile(alpha_MA, beta_MA, discount_factor(iref), ...
                                  yf(iref), sigma_ATM(iref), chi_grid);

    [a0m, a1m, a2m] = local_quad_fit(chi_grid, w_model, CHI_WIN);
    level_model = a0m;
    skew_model  = a1m;
    curv_model  = 2*a2m;

    % Separability self-check: max smile deviation at first/last vs reference.
    w_first = ma_normalized_smile(alpha_MA, beta_MA, discount_factor(1), ...
                                  yf(1), sigma_ATM(1), chi_grid);
    w_last  = ma_normalized_smile(alpha_MA, beta_MA, discount_factor(M), ...
                                  yf(M), sigma_ATM(M), chi_grid);
    invariance_dev = max([max(abs(w_first - w_model)), max(abs(w_last - w_model))]);

    % ===================================================================== %
    %  4) REPORT                                                            %
    % ===================================================================== %
    fprintf('\n=========================================================================\n');
    fprintf('               MA SKEW DIAGNOSTIC (Bachelier implied vol)                \n');
    fprintf('=========================================================================\n');
    fprintf('Calibrated MA:  alpha = %.4f,  beta = %.4f,  beta/alpha = %.4f\n', ...
            a, b, b/a);
    fprintf('Asym. Laplace skewness (distribution, theory): gamma1 = %+.4f\n', gamma1_dist);
    fprintf('  (sign: beta<alpha => positive skew / heavier right tail)\n');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('MA model normalized smile (in chi):  ATM level = %.4f (target ~1.0)\n', level_model);
    fprintf('                                     ATM SKEW  = %+.5f\n', skew_model);
    fprintf('                                     curvature = %+.5f\n', curv_model);
    fprintf('Separability self-check: max |w(t) - w(t_ref)| over chi = %.2e\n', invariance_dev);
    fprintf('  (should be ~0: confirms the MA smile is maturity-independent in chi)\n');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('%-15s | %-10s | %-10s | %-10s | %-12s\n', ...
            'Expiry', 'ATM level', 'ATM skew', 'curvature', 'skew-model');
    fprintf('-------------------------------------------------------------------------\n');
    for i = 1:M
        fprintf('%-15s | %-10.4f | %+-10.5f | %+-10.5f | %+-12.5f\n', ...
            string(expiries(i), 'yyyy-MM-dd'), level_mkt(i), skew_mkt(i), ...
            curv_mkt(i), skew_mkt(i) - skew_model);
    end
    fprintf('-------------------------------------------------------------------------\n');
    valid_sk = ~isnan(skew_mkt);
    fprintf('Market skew (mean +/- std over maturities): %+.5f +/- %.5f\n', ...
            mean(skew_mkt(valid_sk)), std(skew_mkt(valid_sk)));
    fprintf('Model skew: %+.5f   |   mean gap (mkt - model): %+.5f\n', ...
            skew_model, mean(skew_mkt(valid_sk)) - skew_model);
    fprintf('=========================================================================\n\n');

    % ===================================================================== %
    %  5) PLOTS                                                             %
    % ===================================================================== %
    figure('Name', 'MA skew diagnostic');

    % (a) Collapse plot: normalized market smiles vs chi, with the MA curve.
    subplot(1, 2, 1); hold on; grid on;
    cmap = lines(M);
    for i = 1:M
        v = ~isnan(chi_mkt_all(i, :)) & ~isnan(w_mkt_all(i, :));
        if any(v)
            plot(chi_mkt_all(i, v), w_mkt_all(i, v), 'o', ...
                'Color', cmap(i, :), 'MarkerFaceColor', cmap(i, :), ...
                'MarkerSize', 4, 'DisplayName', string(expiries(i), 'MMM-yy'));
        end
    end
    plot(chi_grid, w_model, 'k-', 'LineWidth', 2, 'DisplayName', 'MA model');
    xline(0, ':k', 'HandleVisibility', 'off');
    xlabel('\chi = (K-F)/(\sigma_{ATM}\surd t)');
    ylabel('\sigma_{imp}/\sigma_{ATM}');
    title('Normalized smile collapse (market vs MA)');
    legend('Location', 'best'); xlim([-CHI_WIN*1.5, CHI_WIN*1.5]);

    % (b) Skew term structure: market per maturity vs constant model skew.
    subplot(1, 2, 2); hold on; grid on;
    plot(expiries, skew_mkt, 'o-', 'LineWidth', 1.5, 'DisplayName', 'market');
    yline(skew_model, 'r--', 'LineWidth', 1.5, 'DisplayName', 'MA model (flat)');
    xtickformat('MMM-yy');
    xlabel('Expiry'); ylabel('ATM skew  (dw/d\chi)');
    title('Skew term structure: market vs MA');
    legend('Location', 'best');

    % ----------------------------- output --------------------------------
    skew_report = struct( ...
        'chi_grid',    chi_grid, ...
        'w_model',     w_model, ...
        'skew_model',  skew_model, ...
        'level_model', level_model, ...
        'curv_model',  curv_model, ...
        'skew_mkt',    skew_mkt, ...
        'level_mkt',   level_mkt, ...
        'curv_mkt',    curv_mkt, ...
        'gamma1_dist', gamma1_dist, ...
        'invariance_dev', invariance_dev);
end

% ========================================================================= %
%                          LOCAL HELPER FUNCTIONS                           %
% ========================================================================= %

function w = ma_normalized_smile(alpha, beta, df, yf, sigma_ATM, chi_grid)
% Build the MA call prices on chi_grid for a single maturity, then invert them
% to Bachelier implied vol and normalize by sigma_ATM. Returns w(chi).
    chi_row = chi_grid(:).';                       % 1 x NGRID
    c_mod = price_MA([alpha, beta], df, yf, sigma_ATM, chi_row);  % 1 x NGRID
    c_mod = max(c_mod, 0);                          % clamp FFT-free artefacts
    fmk = -chi_row * sigma_ATM * sqrt(yf);          % F - K
    sig_imp = bachelier_iv(c_mod, df, fmk, yf);
    w = sig_imp / sigma_ATM;
end

function [a0, a1, a2] = local_quad_fit(chi, w, chi_win)
% Local quadratic fit w ~ a0 + a1*chi + a2*chi^2 on the window |chi| <= chi_win.
% a1 is the ATM skew (slope), a0 the ATM level, 2*a2 the curvature.
    chi = chi(:); w = w(:);
    m = abs(chi) <= chi_win & ~isnan(w);
    if nnz(m) < 3
        a0 = NaN; a1 = NaN; a2 = NaN; return;
    end
    X = [ones(nnz(m), 1), chi(m), chi(m).^2];
    coef = X \ w(m);
    a0 = coef(1); a1 = coef(2); a2 = coef(3);
end
