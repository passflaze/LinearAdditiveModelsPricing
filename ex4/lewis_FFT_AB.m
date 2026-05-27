function call_price = lewis_FFT_AB(x_money, t1, t2, k, eta, sigma_t1, sigma_t2, fwd_factor)
% LEWIS_FFT_AB  Bachelier-Lewis price of E[(Delta - x)^+] with
% Delta := f_{t2,t2} - f_{t1,t2} (AB-t2 increment) under Additive Bachelier
% (alpha = 1/2), consistent with Lemma 2 of Forward.pdf.
%
% Under Lemma 2 strict: f_{t1,t2} = fwd_factor * f_{t1,t1} with
% fwd_factor = B(0,t1)/B(0,t2). The CF of the increment of the AB-t2 process
% from s=t1 to s=t2 is
%       phi_Delta(u) = phi_{t2}(u) / phi_{t1}(fwd_factor * u)
% with phi_{t_j} as in Eq.(4) of Baviera-Massaria 2026.
%
% Inputs (column vectors / scalars)
%   x_money    : N x 1 moneyness K2 - F(t1,t2) in $
%   t1, t2     : maturities, with t1 < t2
%   k, eta     : AB parameters (constant in time, Prop.2.1)
%   sigma_t1   : scale at t1, = sigma_ATM(t1)/I_0   (Eq.15)
%   sigma_t2   : scale at t2, = sigma_ATM(t2)/I_0
%   fwd_factor : (optional, default 1) Lemma 2 rescaling B(0,t1)/B(0,t2)
%
% Output
%   call_price : N x 1 undiscounted price E[(Delta - x_money)^+] in $.
%                Multiply by D(t1,t2) externally to get the inner Call at t1.

if nargin < 8 || isempty(fwd_factor)
    fwd_factor = 1;
end

x_money = x_money(:);

% --- FFT grid -----------------------------------------------------------
% Operate on UNNORMALIZED dollar increments (sigma ~ 10-20$), not normalized chi differently from call_ab_FFT.
% We need a large spatial domain to avoid truncation and a very fine dx to 
% resolve the narrow CF peak and the near-pole of the integrand at u=0.
M  = 16;
dz = 0.05;
N  = 2^M;
dx = 2*pi / (N*dz);

z1 = -dz * (N - 1) / 2;
x1 = -dx * (N - 1) / 2;

j  = 0:N-1;
zk = z1 + dz * j;          % output grid (moneyness, $)
xk = x1 + dx * j;          % Fourier integration grid

% --- Lewis contour 'a' inside the (-p+, p-) strip of phi_Delta ----------
% phi_t2(u) strip in u:                |Im u| < p+/(sigma_t2*sqrt(t2)).
% phi_t1(fwd_factor*u) strip in u:     |Im u| < p+/(fwd_factor*sigma_t1*sqrt(t1)).
% Take the binding (tightest) side under Lemma 2.
p_plus    = eta + sqrt(eta^2 + 1/k);
strip_t2  = p_plus / (sigma_t2 * sqrt(t2));
strip_t1e = p_plus / (fwd_factor * sigma_t1 * sqrt(t1));
strip     = min(strip_t2, strip_t1e);
a         = max(-0.5, -0.45 * strip);

% --- Ratio characteristic function for the increment Delta --------------
% Lemma 2: f_{t1,t2} = fwd_factor * f_{t1,t1} -> CF rescaled in the t1 leg.
phi_t1    = charateristic_function_AB(t1, k, eta, sigma_t1);
phi_t2    = charateristic_function_AB(t2, k, eta, sigma_t2);
phi_Delta = @(u) phi_t2(u) ./ phi_t1(fwd_factor * u);

% --- Lewis integrand (Eq.8, a < 0 so R_a = 0) ---------------------------
int = @(csi) phi_Delta(csi + 1i*a) ./ (1i*csi - a).^2;

fj_raw = arrayfun(int, xk);
n = sum(~isfinite(fj_raw));
fj_raw(~isfinite(fj_raw)) = 0; % Prevent NaNs from exploding the FFT at the tails due to 0/0 issues in CF
fj    = fj_raw .* exp(-1i * z1 * dx .* j);
f_hat = exp(a*zk) .* dx .* exp(-1i * x1 * zk) .* fft(fj);

call_price = interp1(zk, f_hat, x_money, 'spline');
call_price = real(call_price / (2*pi));

% No-arb lower bound: E[(Delta-x)+] >= max(0,-x)  (Jensen, E[Delta]=0).
% Clips spline oscillations in the OTM tail; also makes Put via parity
% (Ptilde = Ctilde + x) automatically respect its own bound max(0, x).
call_price = max(call_price, max(0, -x_money));

end
