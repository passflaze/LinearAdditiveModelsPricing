function call_price_modelAB = call_AB_FFT(x, k, eta, I_0)
% CALL_AB_FFT  Normalized AB call G(chi) via single-shift Lewis-FFT.
%   x = moneyness chi at which to interpolate; params [k; eta]; I_0 optional
%   (pass it in during calibration to skip one I0_AB FFT per evaluation).

M=14;
dz=0.005;
N = 2^M;
dx = 2 * pi / (N * dz);

% Lower bounds for frequency (z) and space (x) grids
z1 = -dz * (N - 1) / 2;
x1 = -dx * (N - 1) / 2;

j = 0:N-1;
zk = z1 + dz * j;
xk = x1 + dx * j;

if nargin < 4 || isempty(I_0)
    I_0 = I0_AB(0, [k; eta]);
end

% Lewis contour a in (-p+*I_0, p-*I_0), p+- = +-eta + sqrt(eta^2 + 1/k).
% Default -0.5; pulled in for narrow strips, which also caps exp(a*zk) as k->0.
p_plus = eta + sqrt(eta^2 + 1/k);
a = -min(0.45 * p_plus * I_0, 0.5);

% CF at scale 1/I_0 (so sigma_t = 1/I_0); params = [k; eta].
phi_vals = cf_AB(xk + 1i*a, [k; eta], 1/I_0);
phi_vals = phi_vals(:).';
fj = phi_vals ./ (1i*xk - a).^2 .* exp(-1i * z1 * dx .* j);
f_hat = exp(a*zk) .* dx .* exp(-1i * x1 * zk) .* fft(fj);

call_price_modelAB = interp1(zk, f_hat, x, 'spline');
call_price_modelAB = real(call_price_modelAB/(2*pi));   % f_hat = 2*pi * C_Lewis


end












