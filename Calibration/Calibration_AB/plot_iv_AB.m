function plot_iv_AB(k_AB, eta_AB, discount_factor, yf, sigma_ATM, ...
                    moneyness_modified, c_mkt_calibration, expiries)
% PLOT_IV_AB  Bachelier implied-volatility diagnostic for the Additive
% Bachelier (AB) model: calibrated model vs market.
%
%   Two panels (analogous to check_skew_MA):
%     (a) Normalized smile collapse: w(chi) = sigma_imp/sigma_ATM for every
%         maturity's market quotes (points) against the single AB model curve.
%         By AB separability the normalized model smile is maturity-independent
%         in chi (same property exploited by price_AB / moneyness_generator),
%         so all maturities should collapse onto one curve.
%     (b) Absolute Bachelier IV vs dollar moneyness (F-K) for a representative
%         (median) maturity: market points vs the AB model line.
%
% INPUTS (conventions as in calibrate_surface / check_skew_MA)
%   k_AB, eta_AB        : calibrated AB parameters (scalars), params = [k; eta]
%   discount_factor     : (M x 1) discount factors per maturity
%   yf                  : (M x 1) year fractions per maturity
%   sigma_ATM           : (M x 1) Bachelier ATM implied vols per maturity
%   moneyness_modified  : (M x N) normalized moneyness chi (NaN where no data)
%   c_mkt_calibration   : (M x N) market call prices (NaN where no data)
%   expiries            : (M x 1) datetime expiry dates (for labels)
%
% Implied vols are obtained with the shared Newton inverter Utilities/bachelier_iv.

    params = [k_AB; eta_AB];

    sigma_ATM = sigma_ATM(:); yf = yf(:); discount_factor = discount_factor(:);
    M = numel(yf);

    CHI_WIN = 3.0;       % chi range for the model curve / x-limits
    NGRID   = 401;

    % ---------------------------------------------------------------------
    % 1) MARKET: invert every maturity's quotes to Bachelier IV
    % ---------------------------------------------------------------------
    w_mkt_all   = nan(size(moneyness_modified));   % normalized sigma_imp/sigma_ATM
    sig_mkt_all = nan(size(moneyness_modified));   % absolute sigma_imp
    chi_all     = moneyness_modified;

    for i = 1:M
        chi_i = moneyness_modified(i, :);
        c_i   = c_mkt_calibration(i, :);
        valid = ~isnan(chi_i) & ~isnan(c_i) & (c_i > 0);
        if ~any(valid), continue; end

        chi_v   = chi_i(valid);
        fmk     = -chi_v * sigma_ATM(i) * sqrt(yf(i));    % F - K
        sig_imp = bachelier_iv(c_i(valid), discount_factor(i), fmk, yf(i)).';

        sig_mkt_all(i, valid) = sig_imp;
        w_mkt_all(i, valid)   = sig_imp / sigma_ATM(i);
    end

    % ---------------------------------------------------------------------
    % 2) AB MODEL normalized smile on a dense chi grid (maturity-independent)
    % ---------------------------------------------------------------------
    chi_grid = linspace(-CHI_WIN, CHI_WIN, NGRID);
    iref = max(1, round(M/2));                              % representative maturity

    w_model = ab_normalized_smile(params, discount_factor(iref), yf(iref), ...
                                  sigma_ATM(iref), chi_grid);

    % separability self-check (first vs last maturity vs reference)
    w_first = ab_normalized_smile(params, discount_factor(1), yf(1), sigma_ATM(1), chi_grid);
    w_last  = ab_normalized_smile(params, discount_factor(M), yf(M), sigma_ATM(M), chi_grid);
    invariance_dev = max([max(abs(w_first - w_model)), max(abs(w_last - w_model))], [], 'omitnan');

    % ---------------------------------------------------------------------
    % 3) PLOTS
    % ---------------------------------------------------------------------
    figure('Name', 'AB implied-volatility diagnostic', 'Color', 'w');
    cmap = lines(M);

    % (a) Normalized smile collapse: market points (all maturities) + AB curve
    subplot(1, 2, 1); hold on; grid on;
    for i = 1:M
        v = ~isnan(chi_all(i, :)) & ~isnan(w_mkt_all(i, :));
        if any(v)
            plot(chi_all(i, v), w_mkt_all(i, v), 'o', ...
                'Color', cmap(i, :), 'MarkerFaceColor', cmap(i, :), ...
                'MarkerSize', 4, 'DisplayName', string(expiries(i), 'MMM-yy'));
        end
    end
    plot(chi_grid, w_model, 'k-', 'LineWidth', 2, 'DisplayName', 'AB model');
    xline(0, ':k', 'HandleVisibility', 'off');
    xlabel('\chi = (K-F)/(\sigma_{ATM}\surd t)');
    ylabel('\sigma_{imp}/\sigma_{ATM}');
    title(sprintf('AB normalized IV smile (separability dev = %.1e)', invariance_dev));
    legend('Location', 'best'); xlim([-CHI_WIN, CHI_WIN]);

    % (b) Absolute IV vs dollar moneyness, representative maturity
    subplot(1, 2, 2); hold on; grid on;
    chi_r = chi_all(iref, :);
    sig_r = sig_mkt_all(iref, :);
    vr    = ~isnan(chi_r) & ~isnan(sig_r);
    fmk_r = -chi_r(vr) * sigma_ATM(iref) * sqrt(yf(iref));
    plot(fmk_r, sig_r(vr), 'o', 'Color', cmap(iref, :), ...
        'MarkerFaceColor', cmap(iref, :), 'MarkerSize', 5, 'DisplayName', 'market');

    sig_model = ab_absolute_iv(params, discount_factor(iref), yf(iref), sigma_ATM(iref), chi_grid);
    fmk_grid  = -chi_grid * sigma_ATM(iref) * sqrt(yf(iref));
    plot(fmk_grid, sig_model, 'k-', 'LineWidth', 2, 'DisplayName', 'AB model');
    yline(sigma_ATM(iref), ':k', 'HandleVisibility', 'off');
    xlabel('F - K  ($)'); ylabel('\sigma_{imp}  (Bachelier, $)');
    title(sprintf('AB IV smile @ %s', string(expiries(iref), 'MMM-yyyy')));
    legend('Location', 'best');
end

% ===================================================================== %
%                           LOCAL HELPERS                               %
% ===================================================================== %
function w = ab_normalized_smile(params, df, yf, sigma_ATM, chi_grid)
% AB call prices on chi_grid (single maturity), inverted to Bachelier IV and
% normalized by sigma_ATM. Returns w(chi) as a row aligned with chi_grid.
    sig = ab_absolute_iv(params, df, yf, sigma_ATM, chi_grid);
    w = sig / sigma_ATM;
end

function sig = ab_absolute_iv(params, df, yf, sigma_ATM, chi_grid)
% Absolute Bachelier IV of the AB model on chi_grid for one maturity.
    chi_row = chi_grid(:).';                              % 1 x NGRID
    c_mod = price_AB(params, df, yf, sigma_ATM, chi_row); % 1 x NGRID
    c_mod = max(c_mod, 0);                                % clamp tiny FFT artefacts
    fmk = -chi_row * sigma_ATM * sqrt(yf);                % F - K
    sig = bachelier_iv(c_mod, df, fmk, yf).';             % row, aligned with chi_grid
end
