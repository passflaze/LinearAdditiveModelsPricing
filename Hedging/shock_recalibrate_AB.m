function state = shock_recalibrate_AB(market, params0, dvol)
% SHOCK_RECALIBRATE_AB  Sticky-strike vol shock + AB re-calibration.
%
%   Implements the "model-consistent" vega scenario: every market option keeps
%   its strike and its forward/discount, but its Bachelier implied vol is
%   shifted by a PARALLEL absolute amount dvol (sticky-strike). The shocked
%   surface is re-priced and the AB model is RE-CALIBRATED to it, so the shape
%   parameters (k, eta) are free to move. Contrast with the analytic
%   scale-bump vega (self-similar / sticky-moneyness), where (k, eta) are held
%   fixed by construction.
%
%   Pipeline (per maturity i, on the SAME quoted points used in calibration):
%     1. invert the calibration call prices to Bachelier IV     (bachelier_iv)
%     2. shock:  IV -> IV + dvol                                (sticky-strike)
%     3. re-price the shocked calls in closed form              (Bachelier)
%     4. shock the ATM vol:  sigma_ATM -> sigma_ATM + dvol
%     5. recompute the normalized moneyness chi with the shocked ATM vol
%     6. re-calibrate AB via fmincon, WARM-STARTED at params0   (low jitter)
%
%   The valid-point mask is derived from the (dvol-independent) base IV
%   inversion, so the +dvol and -dvol calibrations fit the SAME target set:
%   this keeps the central finite difference clean.
%
% INPUTS:
%   market  - (struct) from run_ex2: needs .forward .strikes .sigma_ATM .yf
%                      .discount_factor .moneyness_modified .c_mkt_calibration
%   params0 - (2x1) base AB params [k; eta] (warm-start for the optimizer)
%   dvol    - (scalar) signed sticky-strike vol shift (absolute Bachelier vol)
%
% OUTPUT:
%   state   - (struct) .params_AB       recalibrated [k; eta]
%                      .sigma_ATM        shocked ATM vol vector (M x 1)
%                      .fval .exitflag   calibration diagnostics

    forward = market.forward(:);          % M x 1
    strikes = market.strikes(:).';        % 1 x N
    yf      = market.yf(:);               % M x 1
    df      = market.discount_factor(:);  % M x 1
    c_base  = market.c_mkt_calibration;   % M x N (NaN off-grid)
    [M, N]  = size(c_base);

    % --- Broadcast helpers (M x N) ---------------------------------------
    FmK    = forward - strikes;           % F - K
    B_mat  = repmat(df, 1, N);
    t_mat  = repmat(yf, 1, N);

    % --- 1. Base Bachelier IV on the quoted points -----------------------
    IV_base = reshape(bachelier_iv(c_base, B_mat, FmK, t_mat), M, N);

    % Valid set: quoted AND invertible. dvol-independent -> shared by +/-dvol.
    valid = ~isnan(c_base) & isfinite(IV_base);

    % --- 2-3. Sticky-strike shock + closed-form Bachelier re-price -------
    IV_shk = IV_base + dvol;
    v      = IV_shk .* sqrt(t_mat);
    d      = FmK ./ v;
    Phi    = 0.5 * erfc(-d / sqrt(2));
    phi    = exp(-0.5 * d.^2) / sqrt(2*pi);
    c_shk  = B_mat .* (FmK .* Phi + v .* phi);

    c_shk(~valid) = NaN;                  % keep only the shared target set

    % --- 4. Sticky-strike also lifts the ATM vol -------------------------
    sigma_ATM_shk = market.sigma_ATM(:) + dvol;   % M x 1

    % --- 5. Normalized moneyness with the shocked ATM vol ----------------
    chi_full = (strikes - forward) ./ (sigma_ATM_shk .* sqrt(yf));   % M x N
    chi_shk  = NaN(M, N);
    chi_shk(valid) = chi_full(valid);

    % --- 6. Re-calibrate AB, warm-started at params0 ---------------------
    options_AB = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point', ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance',       1e-10, ...
        'MaxFunctionEvaluations', 5000);
    lb_AB = [1e-3, -1.5];
    ub_AB = [5.0,   1.5];
    x0    = params0(:).';                 % warm start (row, as fmincon expects)

    obj = @(x) objective_function_AB(x, df, yf, sigma_ATM_shk, chi_shk, c_shk);
    [x_opt, fval, exitflag] = fmincon(obj, x0, [], [], [], [], lb_AB, ub_AB, [], options_AB);

    state = struct('params_AB', [x_opt(1); x_opt(2)], ...
                   'sigma_ATM',  sigma_ATM_shk, ...
                   'fval', fval, 'exitflag', exitflag);
end
