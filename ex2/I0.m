function I = I0(x,k, eta)
% Computes call prices using FFT method
% Inputs/output to be completed

M=14;
dz=0.005;
N = 2^M;
dx = 2 * pi / (N * dz);

% Lower bounds for frequency (z) and space (x) grids
z1 = -dz * (N - 1) / 2;
x1 = -dx * (N - 1) / 2;
a=-1/2;
j = 0:N-1;
zk = z1 + dz * j;
xk = x1 + dx * j;


%to be completed