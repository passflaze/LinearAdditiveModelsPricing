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

a=-1/2; % countour offset for the Fourier inversion, in the strip of regularity of the cf

phi=charateristic_function_AB(1,k,eta,1);

int=@(csi) phi(csi+1i*a)./(1i*csi-a).^2;

fj= arrayfun(int,xk) .* exp(-1i * z1 * dx .*j);

f_hat= exp(a*zk) .* dx .* exp(-1i * x1 * zk) .* fft(fj);

I= interp1(zk, f_hat, x, 'spline');

I= real(I/sqrt(2*pi));

end