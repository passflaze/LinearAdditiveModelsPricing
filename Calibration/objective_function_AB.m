function sse = objective_function_AB(x, df, yf, sigma_ATM, mon_mod, c_mkt)
% OBJECTIVE_FUNCTION_AB Computes scalar SSE between AB model and market prices.
%
% Robust wrapper around the Additive Bachelier pricing engine (price_AB) for
% calibration via fmincon. Mirrors objective_function_MA / _GL: it prices the
% whole (M x N) surface in one shot, penalizes model failures, and sums the
% squared residuals over the quoted points only. NaNs in c_mkt (strikes
% outside the band or without market data) are ignored.
%
% INPUTS:
%   x         : 2-element parameter vector [k, eta]
%   df        : (M x 1) discount factors
%   yf        : (M x 1) year fractions (time to maturity)
%   sigma_ATM : (M x 1) ATM Bachelier implied vols
%   mon_mod   : (M x N) modified market moneyness (chi), NaNs allowed
%   c_mkt     : (M x N) market call calibration prices, NaNs allowed
%
% OUTPUT:
%   sse       : scalar SSE to be minimized

    k   = x(1);
    eta = x(2);

    % --- Step 1: Execute Pricing Engine ---
    c_mod = price_AB(x, df, yf, sigma_ATM, mon_mod);

    % --- Step 2: STRICT PENALTY CHECK ---
    % Se c_mod è NaN (fallimento del modello) MA c_mkt è un numero valido, scatta la sanzione.
    if any(isnan(c_mod(:)) & ~isnan(c_mkt(:)))
        fprintf('[PENALTY] Il modello ha restituito NaN su un target valido di mercato. k: %.4f, eta: %.4f\n', k, eta);
        sse = 1e10 + 1e5 * ((k - 1)^2 + eta^2);
        return;
    end

    % --- Step 3: SAFE SSE CALCULATION ---
    % Sottraiamo le matrici. Se in una cella c'è un NaN (es. l'opzione di mercato non c'era),
    % la differenza rimarrà semplicemente NaN.
    diff = c_mkt - c_mod;

    % Creiamo una maschera che seleziona SOLO i punti in cui la differenza è un numero vero.
    valid_idx = ~isnan(diff);

    % Controllo anti-collasso (evita il finto "SSE = 0" che inganna fmincon)
    if sum(valid_idx(:)) == 0 || all(c_mod(valid_idx) == 0)
        fprintf('[PENALTY] La superficie del modello è collassata a zero. k: %.4f, eta: %.4f\n', k, eta);
        sse = 1e10 + 1e5 * ((k - 1)^2 + eta^2);
        return;
    end

    % Estraiamo solo i residui validi (diventa un vettore) ed eleviamo al quadrato
    residuals = diff(valid_idx);
    sse = sum(residuals.^2);

end
