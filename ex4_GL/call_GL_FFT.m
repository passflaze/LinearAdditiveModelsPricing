function call_price = call_GL_FFT(x_money, t1, t2, alpha, beta, sigma_t1, sigma_t2)
% CALL_GL_FFT  Bachelier-Lewis FFT price of E[(Delta - x)^+] with
% Delta := f_{t2} - f_{t1} under the Generalized Logistic additive model.
%
% GL is additive but NOT Levy, so Delta is NOT GL-distributed: its CF is
% the ratio   phi_Delta(u) = phi_{t2}(u) / phi_{t1}(u)
% with phi_{t_j}(u) = cf_GL( u * sigma_{t_j} * sqrt(t_j), alpha, beta ).
% The Bachelier-Lewis formula is applied directly to phi_Delta on a
% contour 'a' inside the analyticity strip: in u-space the strip is
%   ( -alpha/(sigma_t2 sqrt(t2)),  beta/(sigma_t2 sqrt(t2)) )
% (t2 is the narrower one, since sigma_t*sqrt(t) is increasing in t).
%
% Inputs
%   x_money    : N x 1 moneyness K2 - F(t1,t2) in $
%   t1, t2     : maturities, with t1 < t2
%   alpha, beta: GL shape parameters
%   sigma_t1   : scale at t1, = sigma_ATM(t1)/I_0
%   sigma_t2   : scale at t2, = sigma_ATM(t2)/I_0
%
% Output
%   call_price : N x 1 UNDISCOUNTED price E[(Delta - x_money)^+] in $.
%                Multiply by B(t1,t2) externally to get the inner Call at t1.

x_money = x_money(:);

% --- FFT grid (unnormalized dollar increments) --------------------------
M  = 16;
dz = 0.005;
N  = 2^M;
dx = 2*pi / (N*dz);

z1 = -dz * (N - 1) / 2;
x1 = -dx * (N - 1) / 2;

j  = 0:N-1;
zk = z1 + dz * j;          % output grid (moneyness, $)
xk = x1 + dx * j;          % Fourier integration grid

% --- Lewis contour 'a' inside the (-alpha/.., beta/..) strip ------------
% Intersection of the two strips is dominated by t2 (narrower).
strip_lo = -alpha / (sigma_t2 * sqrt(t2));   % < 0
a        = max(-0.5, 0.45 * strip_lo);       % a < 0  -> R_a = 0

% --- Ratio CF for the increment Delta -----------------------------------
phi_t1    = @(u) cf_GL(u .* (sigma_t1 * sqrt(t1)), alpha, beta);
phi_t2    = @(u) cf_GL(u .* (sigma_t2 * sqrt(t2)), alpha, beta);
phi_Delta = @(u) phi_t2(u) ./ phi_t1(u);

% --- Lewis integrand (Baviera-Manzoni eq.8, a<0 -> R_a=0) ---------------
%   integrand(xi) = phi_Delta(xi + i a) / (i xi - a)^2
int = @(csi) phi_Delta(csi + 1i*a) ./ (1i*csi - a).^2;

fj_raw = arrayfun(int, xk);
fj_raw(~isfinite(fj_raw)) = 0;               % kill 0/0 at the tails
fj    = fj_raw .* exp(-1i * z1 * dx .* j);
f_hat = exp(a*zk) .* dx .* exp(-1i * x1 * zk) .* fft(fj);

call_price = interp1(zk, f_hat, x_money, 'spline');
call_price = real(call_price / (2*pi));

% No-arb floor: E[(Delta - x)^+] >= max(0, -x)  (Jensen, E[Delta]=0).
call_price = max(call_price, max(0, -x_money));

end