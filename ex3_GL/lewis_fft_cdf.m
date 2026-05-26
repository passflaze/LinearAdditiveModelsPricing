function [CDF_grid, zk] = lewis_fft_cdf(cf_increment, alpha, beta, ...
                                       sigma_s, s, sigma_t, t)

    M  = 16;
    N  = 2^M;
    dz = 0.005;
    dx = 2*pi / (N * dz);
    z1 = -dz * (N - 1) / 2;
    x1 = -dx * (N - 1) / 2;
    j  = 0:N-1;
    zk = z1 + dz * j;
    xk = x1 + dx * j;

    p_minus_t = alpha / (sigma_t * sqrt(t));
    p_plus_t  = beta  / (sigma_t * sqrt(t));
    a_neg = -p_minus_t / 2;
    a_pos =  p_plus_t  / 2;

    int_fn = @(u, a) cf_increment(u - 1i*a, alpha, beta, ...
                                  sigma_s, s, sigma_t, t) ./ (1i*u + a);

    fj_neg    = int_fn(xk, a_neg) .* exp(-1i * z1 * dx .* j);
    I_hat_neg = real( dx .* exp(-1i * x1 * zk) .* fft(fj_neg) );

    fj_pos    = int_fn(xk, a_pos) .* exp(-1i * z1 * dx .* j);
    I_hat_pos = real( dx .* exp(-1i * x1 * zk) .* fft(fj_pos) );

    CDF_grid = zeros(size(zk));

    idx_sx = zk < -1;
    CDF_grid(idx_sx) = - (exp(-a_neg * zk(idx_sx)) / (2*pi)) .* I_hat_neg(idx_sx);

    idx_dx = zk > 1;
    CDF_grid(idx_dx) = 1 - (exp(-a_pos * zk(idx_dx)) / (2*pi)) .* I_hat_pos(idx_dx);

    idx_blend = (zk >= -1) & (zk <= 1);
    z_b  = zk(idx_blend);
    c_sx = - (exp(-a_neg * z_b) / (2*pi)) .* I_hat_neg(idx_blend);
    c_dx = 1 - (exp(-a_pos * z_b) / (2*pi)) .* I_hat_pos(idx_blend);
    w    = 1 ./ (1 + exp(-5 * z_b));
    CDF_grid(idx_blend) = c_sx .* (1 - w) + c_dx .* w;

end