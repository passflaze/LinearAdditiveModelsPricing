function I = I0(x,k, eta)
% Computes I0 using FFT method
% Inputs/output to be completed

dz=0.005;
M=14;
N = 2^M;
dx = 2 * pi / (N * dz);

% Lower bounds for frequency (z) and space (x) grids
z1 = -dz * (N - 1) / 2;
x1 = -dx * (N - 1) / 2;
j = 0:N-1;
zk = z1 + dz * j;
xk = x1 + dx * j;

% Lewis contour: must lie inside (-p+, p-) with p+- = +-eta + sqrt(eta^2 + 1/k)
% (sigma_t = 1 inside the CF here, so the f_t-strip coincides with the zeta-strip).
% Default -1/2 if safely inside; otherwise pull it back to a fraction of p+ with margin.
p_plus = eta + sqrt(eta^2 + 1/k);
a = -min(0.45 * p_plus, 0.5);

phi=charateristic_function_AB(1,k,eta,1);

% vectorised CF call: phi(...) returns column N x 1; reshape to row for FFT
phi_vals = phi(xk + 1i*a);
phi_vals = phi_vals(:).';
fj = phi_vals ./ (1i*xk - a).^2 .* exp(-1i * z1 * dx .* j);

f_hat= exp(a*zk) .* dx .* exp(-1i * x1 * zk) .* fft(fj);

I= interp1(zk, f_hat, x, 'spline');

I= real(I/sqrt(2*pi));

end