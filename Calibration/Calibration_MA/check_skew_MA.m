function skew_report = check_skew_MA(alpha_MA, beta_MA, discount_factor, yf, ...
                                     sigma_ATM, moneyness_modified, ...
                                     c_mkt_calibration, expiries)
% CHECK_SKEW_MA  Post-calibration skew diagnostic for the Minimal Additive model.
%
%   Verifies a calibrated MA fit using the Bachelier implied-volatility SKEW
%   as the primary check. This is the natural diagnostic for MA because:
%
%     1) The MA driver zeta_t is an asymmetric Laplace with rates (alpha,beta)
%        (cf. price_MA.m: density ~ e^{alpha x} for x<0, ~ e^{-beta x} for x>0).
%        Its skewness is SCALE-INVARIANT and depends ONLY on the ratio
%        beta/alpha -- which is exactly the single identifiable parameter of
%        the model (the (alpha,beta)->(s*alpha,s*beta) gauge invariance
%        documented in run_exercise_2.m). So checking the skew == checking the
%        one number you actually calibrated.
%
%     2) Under the additive-Bachelier separability of the implied vol
%        (Baviera & Massaria 2026, the property used in moneyness_generator.m),
%        the MA implied-vol smile written in normalized moneyness
%        chi = (K-F)/(sigma_ATM*sqrt(t)) is MATURITY-INDEPENDENT. Hence the
%        per-maturity market smiles, once normalized, should collapse onto a
%        single curve whose ATM slope equals the model skew. Strong, falsifiable
%        prediction: if market skews drift with maturity, a single-beta MA
%        cannot match them (model limitation, not a calibration bug).
%
%   SKEW METRIC. For each maturity we invert market call prices to the
%   Bachelier implied vol sigma_imp(chi), form the dimensionless normalized
%   smile w(chi) = sigma_imp(chi)/sigma_ATM, and fit locally (|chi| <= CHI_WIN)
%       w(chi) ~ a0 + a1*chi + a2*chi^2.
%   The ATM SKEW is a1 (= dw/dchi at chi=0); a0 is the ATM level (~1, sanity
%   check on the cascade) and 2*a2 the curvature. The MA model skew is measured
%   identically on a dense chi grid, so the comparison is apples-to-apples.
%
%   INPUTS (conventions as in run_exercise_2.m)
%     alpha_MA, beta_MA   : calibrated MA rate parameters (scalars)
%     discount_factor     : (M x 1) discount factors per maturity
%     yf                  : (M x 1) year fractions per maturity
%     sigma_ATM           : (M x 1) Bachelier ATM implied vols
%     moneyness_modified  : (M x N) normalized moneyness chi (NaN where no data)
%     c_mkt_calibration   : (M x N) target market call prices (NaN where no data)
%     expiries            : (M x 1) datetime expiry dates (for labels)
%
%   OUTPUT
%     skew_report : struct with fields
%        .chi_grid        dense chi grid used for the model curve
%        .w_model         model normalized smile w(chi) on chi_grid
%        .skew_model      MA model ATM skew (a1)
%        .level_model     MA model ATM level (a0, should be ~1)
%        .curv_model      MA model ATM curvature (2*a2)
%        .skew_mkt        (M x 1) market ATM skew per maturity
%        .level_mkt       (M x 1) market ATM level per maturity
%        .curv_mkt        (M x 1) market ATM curvature per maturity
%        .gamma1_dist     analytic skewness of the asymmetric Laplace (theory)
%
%   Reference: Baviera & Massaria (2026), "The additive Bachelier model",
%   J. Comput. Appl. Math. 487, 117741 (Bachelier IV, Eq. (15); separability).

    % --------------------------- settings --------------------------------
    CHI_WIN = 2.0;     % half-width of the ATM window (in chi) for the local fit
    NGRID   = 401;     % points on the dense model chi grid

    alpha_MA = alpha_MA(:); beta_MA = beta_MA(:);
    sigma_ATM = sigma_ATM(:); yf = yf(:); discount_factor = discount_factor(:);
    M = numel(yf);

    % ===================================================================== %
    %  1) ANALYTIC ASYMMETRIC-LAPLACE SKEWNESS (theory, distribution-level)  %
    % ===================================================================== %
    % For density ~ e^{alpha x} (x<0), ~ e^{-beta x} (x>0):
    %   gamma1 = 2*(1/beta^3 - 1/alpha^3) / (1/alpha^2 + 1/beta^2)^(3/2).
    % Scale-invariant -> function of beta/alpha only. Sign: beta<alpha => right
    % tail heavier => positive skew.
    a = alpha_MA; b = beta_MA;
    gamma1_dist = 2*(1/b^3 - 1/a^3) / (1/a^2 + 1/b^2)^(3/2);

    % ===================================================================== %
    %  2) MARKET SKEW PER MATURITY (invert prices -> normalized IV -> fit)   %
    % ===================================================================== %
    skew_mkt  = nan(M, 1);
    level_mkt = nan(M, 1);
    curv_mkt  = nan(M, 1);

    % store normalized smiles for the collapse plot
    w_mkt_all   = nan(size(moneyness_modified));
    chi_mkt_all = moneyness_modified;

    for i = 1:M
        chi_i = moneyness_modified(i, :);
        c_i   = c_mkt_calibration(i, :);
        valid = ~isnan(chi_i) & ~isnan(c_i) & (c_i > 0);
        if nnz(valid) < 3, continue; end

        chi_v = chi_i(valid);
        c_v   = c_i(valid);

        % Dollar moneyness F-K = -chi * sigma_ATM * sqrt(t).
        fmk = -chi_v * sigma_ATM(i) * sqrt(yf(i));

        % Invert to Bachelier implied vol, then normalize by sigma_ATM.
        sig_imp = bachelier_iv(c_v, discount_factor(i), fmk, yf(i));
        w_v = sig_imp / sigma_ATM(i);

        w_mkt_all(i, valid) = w_v;

        % Local quadratic fit on the ATM window.
        [a0, a1, a2] = local_quad_fit(chi_v, w_v, CHI_WIN);
        level_mkt(i) = a0;
        skew_mkt(i)  = a1;
        curv_mkt(i)  = 2*a2;
    end

    % ===================================================================== %
    %  3) MODEL SKEW (MA smile on a dense chi grid, maturity-independent)    %
    % ===================================================================== %
    chi_grid = linspace(-CHI_WIN*1.5, CHI_WIN*1.5, NGRID);

    % Use a representative (median) maturity; by separability the normalized
    % MA smile is the same for any maturity -- we exploit that below as a check.
    iref = max(1, round(M/2));

    w_model = ma_normalized_smile(alpha_MA, beta_MA, discount_factor(iref), ...
                                  yf(iref), sigma_ATM(iref), chi_grid);

    [a0m, a1m, a2m] = local_quad_fit(chi_grid, w_model, CHI_WIN);
    level_model = a0m;
    skew_model  = a1m;
    curv_model  = 2*a2m;

    % Maturity-invariance self-check: recompute at the first and last maturity
    % and report the max deviation of the normalized smile from the reference.
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
