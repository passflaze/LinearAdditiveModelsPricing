function plot_iv_AB(k_AB, eta_AB, discount_factor, yf, sigma_ATM, ...
                    moneyness_modified, c_mkt_calibration, expiries)
% PLOT_IV_AB  Bachelier IV diagnostic for the AB model: calibrated vs market.
%   (a) Normalized smile collapse w(chi)=sigma_imp/sigma_ATM (market points vs
%       the single, maturity-independent AB curve).
%   (b) Absolute IV vs dollar moneyness (F-K) at a representative maturity.
%   IVs from Utilities/bachelier_iv.
%
% INPUTS:
%   k_AB, eta_AB        : calibrated AB parameters (params = [k; eta])
%   discount_factor, yf, sigma_ATM : (M x 1) per maturity
%   moneyness_modified  : (M x N) chi (NaN where no data)
%   c_mkt_calibration   : (M x N) market call prices (NaN where no data)
%   expiries            : (M x 1) datetime expiry dates (labels)
% OUTPUT:
%   none (draws the normalized-smile collapse and absolute-IV figures)

    params = [k_AB; eta_AB];

    sigma_ATM = sigma_ATM(:); yf = yf(:); discount_factor = discount_factor(:);
    M = numel(yf);

    CHI_WIN = 3.0;       % chi range for the model curve / x-limits
    NGRID   = 401;

    % ---------------------------------------------------------------------
    % 1) MARKET: invert the WHOLE surface to Bachelier IV in one vectorized call
    % ---------------------------------------------------------------------
    N       = size(moneyness_modified, 2);
    chi_all = moneyness_modified;
    sT      = sigma_ATM .* sqrt(yf);                       % (M x 1)
    fmk_all = -moneyness_modified .* sT;                   % (M x N) F - K
    B_all   = repmat(discount_factor, 1, N);              % (M x N)
    t_all   = repmat(yf, 1, N);                            % (M x N)

    valid = ~isnan(moneyness_modified) & ~isnan(c_mkt_calibration) & (c_mkt_calibration > 0);

    sig_mkt_all = nan(size(moneyness_modified));           % absolute sigma_imp
    sig_mkt_all(valid) = bachelier_iv(c_mkt_calibration(valid), B_all(valid), ...
                                      fmk_all(valid), t_all(valid));
    w_mkt_all = sig_mkt_all ./ sigma_ATM;                  % (M x N) normalized (per-row)

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
