function call_price_modelAB = call_AB_FFT(x,k,eta)


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

a=-1/2; % countour offset for the Fourier inversion, in the strip of regularity of the cf

I_0 = I0(0,k,eta);

phi = charateristic_function_AB(1,k,eta,1/I_0); 

int=@(csi) phi(csi+1i*a)./(1i*csi-a).^2;

fj= arrayfun(int,xk) .* exp(-1i * z1 * dx .*j);
% Lewis prefactor e^{x a}: function of the OUTPUT moneyness grid zk
f_hat= exp(a*zk) .* dx .* exp(-1i * x1 * zk) .* fft(fj);

call_price_modelAB = interp1(zk, f_hat, x, 'spline');
% f_hat = 2 pi * C_Lewis  ->  divide by 2 pi to recover the model price G_hat
call_price_modelAB = real(call_price_modelAB/(2*pi));


end












