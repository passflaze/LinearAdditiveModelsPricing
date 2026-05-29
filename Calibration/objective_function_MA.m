function sse = objective_function_MA(x, df, yf, sigma_ATM, mon_mod, c_mkt)
% OBJECTIVE_FUNCTION_MA Computes scalar SSE between model and market prices.
% NaNs in c_mkt are treated as zero contribution to the sum (ignored).
    alpha = x(1);
    beta  = x(2);
    % Generate the candidate option pricing surface matrix
    c_mod = price_MA(x, df, yf, sigma_ATM, mon_mod);

    % --- Step 2: STRICT PENALTY CHECK ---
    % Se c_mod è NaN (fallimento del modello) MA c_mkt è un numero valido, scatta la sanzione.
    if any(isnan(c_mod(:)) & ~isnan(c_mkt(:)))
        fprintf('[PENALTY] Il modello ha restituito NaN su un target valido di mercato. Alpha: %.4f, Beta: %.4f\n', alpha, beta);
        sse = 1e10 + 1e5 * ((alpha - 1)^2 + (beta - 1)^2);
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
        fprintf('[PENALTY] La superficie del modello è collassata a zero. Alpha: %.4f, Beta: %.4f\n', alpha, beta);
        sse = 1e10 + 1e5 * ((alpha - 1)^2 + (beta - 1)^2);
        return;
    end
    
    % Estraiamo solo i residui validi (diventa un vettore) ed eleviamo al quadrato
    residuals = diff(valid_idx);

    sse = sum(residuals.^2);

end