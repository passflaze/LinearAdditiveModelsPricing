function cf_inc = cf_increment_GL(u, alpha, beta, sigma_s, s, sigma_t, t)
% Compute the CF between any pair of reset dates t>s>t0 using the additivity
    
    u_t = u .* (sigma_t * sqrt(t));
    u_s = u .* (sigma_s * sqrt(s));
    
    cf_t = cf_GL(u_t, alpha, beta);
    cf_s = cf_GL(u_s, alpha, beta);
    
    cf_inc = cf_t ./ cf_s;
    
    % --- QUANT PATCH per stabilità numerica ---
    nan_mask = isnan(cf_inc) | isinf(cf_inc);
    cf_inc(nan_mask) = 0;
end