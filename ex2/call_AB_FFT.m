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
    I_0 = I0(0,k,eta);
end

% Lewis contour: sigma_t = 1/I_0 inside the CF, so the f_t-strip is
% (-p+ * I_0, p- * I_0) with p+- = +-eta + sqrt(eta^2 + 1/k).
% Default -1/2 if safely inside; otherwise pull a back with margin.
p_plus = eta + sqrt(eta^2 + 1/k);
a = max(-0.5, -0.45 * p_plus * I_0);

phi = charateristic_function_AB(1,k,eta,1/I_0);

% vectorised CF call: phi(...) returns column N x 1; reshape to row for FFT
phi_vals = phi(xk + 1i*a);
phi_vals = phi_vals(:).';
fj = phi_vals ./ (1i*xk - a).^2 .* exp(-1i * z1 * dx .* j);
% Lewis prefactor e^{x a}: function of the OUTPUT moneyness grid zk
f_hat= exp(a*zk) .* dx .* exp(-1i * x1 * zk) .* fft(fj);

call_price_modelAB = interp1(zk, f_hat, x, 'spline');
% f_hat = 2 pi * C_Lewis  ->  divide by 2 pi to recover the model price G_hat
call_price_modelAB = real(call_price_modelAB/(2*pi));


end












