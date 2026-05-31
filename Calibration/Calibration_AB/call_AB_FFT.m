function call_price_modelAB = call_AB_FFT(x,k,eta,I_0)
% Optional I_0: when calibrating, the caller computes I_0 once per (k,eta)
% and passes it in to skip an extra 2^14 FFT per objective evaluation.


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

% Lewis contour: sigma_t = 1/I_0 inside the CF, so the f_t-strip is
% (-p+ * I_0, p- * I_0) with p+- = +-eta + sqrt(eta^2 + 1/k).
% Keep the stable default 0.5 when the strip is wide; only pull the contour
% in (towards the lower branch point) when the strip is narrow. The cap also
% prevents exp(a*zk) overflow as k -> 0, where p_plus ~ 1/sqrt(k) blows up.
p_plus = eta + sqrt(eta^2 + 1/k);
a = -min(0.45 * p_plus * I_0, 0.5);

% cf_AB with scale_factor = (1/I_0)*sqrt(1) = 1/I_0; params = [k; eta].
phi_vals = cf_AB(xk + 1i*a, [k; eta], 1/I_0);
phi_vals = phi_vals(:).';
fj = phi_vals ./ (1i*xk - a).^2 .* exp(-1i * z1 * dx .* j);
% Lewis prefactor e^{x a}: function of the OUTPUT moneyness grid zk
f_hat= exp(a*zk) .* dx .* exp(-1i * x1 * zk) .* fft(fj);

call_price_modelAB = interp1(zk, f_hat, x, 'spline');
% f_hat = 2 pi * C_Lewis  ->  divide by 2 pi to recover the model price G_hat
call_price_modelAB = real(call_price_modelAB/(2*pi));


end












