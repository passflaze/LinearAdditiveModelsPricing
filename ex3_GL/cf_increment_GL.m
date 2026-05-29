function cf_inc = cf_increment_GL(u, alpha, beta, sigma_s, s, sigma_t, t)
% Compute the CF between any pair of reset dates t>s>t0 using the additivity.
% Caso speciale: se sigma_s = 0 e s = 0, ritorna direttamente cf_GL(u_t).
    u_t = u .* (sigma_t * sqrt(t));

    if sigma_s == 0 && s == 0
        cf_inc = cf_GL(u_t, alpha, beta);
        return;
    end

    u_s = u .* (sigma_s * sqrt(s));
    cf_t = cf_GL(u_t, alpha, beta);
    cf_s = cf_GL(u_s, alpha, beta);
    cf_inc = cf_t ./ cf_s;
% --- QUANT PATCH per stabilità numerica ---
    nan_mask = isnan(cf_inc) | isinf(cf_inc);
    n_bad = nnz(nan_mask);
    if n_bad > 0
        n_nan = nnz(isnan(cf_inc));
        n_inf = nnz(isinf(cf_inc));
        warning('cf_increment_GL:InstabilitaNumerica', ...
            ['Rilevati %d valori non finiti nella CF incrementale ' ...
             '(NaN: %d, Inf: %d) su %d totali. Impostati a 0.'], ...
            n_bad, n_nan, n_inf, numel(cf_inc));
    end
    cf_inc(nan_mask) = 0;
end