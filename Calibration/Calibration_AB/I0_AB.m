function I = I0_AB(x, params)
% I0_AB  Normalization constant I0 = sqrt(2*pi)*E[z^+] of the AB model, via a
%        Lewis-FFT of the (alpha = 0.5) CF at unit scale. Ties market ATM vol to
%        the internal scale: sigma_t = sigma_ATM / I0.
%
% INPUTS
%   x      : evaluation point(s) of I0(.) (use 0 for I0).
%   params : [k; eta]  (k > 0 tempering parameter, eta drift)
% OUTPUT
%   I      : I0 at x (scalar for scalar x).
%
% I0(0) is memoized on (k, eta): callers hit it repeatedly with the same
% parameters and the FFT is deterministic, so a one-slot cache is exact.

    persistent cache_keta cache_I0

    k   = params(1);
    eta = params(2);

    if isscalar(x) && x == 0 && ~isempty(cache_keta) && ...
            cache_keta(1) == k && cache_keta(2) == eta
        I = cache_I0;
        return;
    end

    dz = 0.005;
    M  = 14;
    N  = 2^M;
    dx = 2 * pi / (N * dz);

    % Lower bounds for frequency (z) and space (x) grids
    z1 = -dz * (N - 1) / 2;
    x1 = -dx * (N - 1) / 2;
    j  = 0:N-1;
    zk = z1 + dz * j;
    xk = x1 + dx * j;

    % Lewis contour: must lie inside (-p+, p-) with p+- = +-eta + sqrt(eta^2 + 1/k).
    % Default -1/2 if safely inside; otherwise pull it back to a fraction of p+.
    p_plus = eta + sqrt(eta^2 + 1/k);
    a = -min(0.45 * p_plus, 0.5);

    % CF evaluated at sigma_t = 1, t = 1 (scale_factor = 1); params = [k; eta].
    phi_vals = cf_AB(xk + 1i*a, params, 1);
    phi_vals = phi_vals(:).';
    fj = phi_vals ./ (1i*xk - a).^2 .* exp(-1i * z1 * dx .* j);

    f_hat = exp(a*zk) .* dx .* exp(-1i * x1 * zk) .* fft(fj);

    I = interp1(zk, f_hat, x, 'spline');
    I = real(I/sqrt(2*pi));

    if isscalar(x) && x == 0
        cache_keta = [k, eta];
        cache_I0   = I;
    end
end
