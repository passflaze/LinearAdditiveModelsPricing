function cdf = lewis_fft_cdf(alpha, beta, T1, T2, sigma_T1, sigma_T2, x, flag)
% LEWIS_FFT_CDF  Conditional CDF F_{T2|T1}(x) of the GL log-price increment
% between reset dates T1 and T2, via Lewis-FFT (Baviera-Manzoni eq.13-15).
%
% Two-shift reconstruction: CDF computed twice -- once with a contour
% shift on the negative side of the analyticity strip (accurate on the
% LEFT tail) and once on the positive side (accurate on RIGHT tail) --
% then glued at x = 0.
%
% Inputs
%   alpha, beta : GL shape parameters
%   T1, T2      : reset dates (T1 may be 0 for marginal CDF of f_{T2})
%   sigma_T1    : scale at T1 (= sigma_ATM(T1)/I_0); pass 0 if T1==0
%   sigma_T2    : scale at T2 (= sigma_ATM(T2)/I_0)
%   x           : column vector of log-price increments where to evaluate
%   flag        : 1 -> coarse dx for exotic pricing (unnormalized $)
%                 0 -> fine dx for diagnostics
%
% Output
%   cdf : column vector of CDF values at x
 
phi_T2 = @(u) cf_GL(u .* (sigma_T2 * sqrt(T2)), alpha, beta);
 
if T1 == 0
    phi_cond = phi_T2;
else
    phi_T1   = @(u) cf_GL(u .* (sigma_T1 * sqrt(T1)), alpha, beta);
    phi_cond = @(u) phi_T2(u) ./ phi_T1(u);
end
 
% FFT grid
M = 16;
if flag == 1
    dz = 0.05;     % unnormalized dollar increments
else
    dz = 0.0025;
end
N  = 2^M;
dx = 2*pi / (N*dz);
z1 = -dz * (N-1) / 2;
x1 = -dx * (N-1) / 2;
j  = 0:N-1;
zk = z1 + dz * j;
xk = x1 + dx * j;
 
grid = struct('zk', zk, 'xk', xk, 'x1', x1, 'z1', z1, 'dx', dx, 'j', j);
 
% Analyticity strip in u-space: u in ( -alpha/sT2, beta/sT2 )
sT2   = sigma_T2 * sqrt(T2);
a_neg = -0.49 * alpha / sT2;   % left edge:  good for RIGHT tail, R_a = 1
a_pos =  0.49 * beta  / sT2;   % right edge: good for LEFT tail,  R_a = 0
 
x = x(:);
cdf_right = one_shift(phi_cond, x, a_neg, 1, grid);   % accurate for x > 0
cdf_left  = one_shift(phi_cond, x, a_pos, 0, grid);   % accurate for x < 0
 
% Glue at x = 0
cdf            = cdf_right;
left_mask      = x < 0;
cdf(left_mask) = cdf_left(left_mask);
 
end
 
% ----------------------------------------------------------------------
function P = one_shift(phi_cond, x, a, Ra, g)
% Single FFT reconstruction of the CDF (Baviera-Manzoni eq.13-15):
%   P(x) = R_a - (e^{a x}/(2 pi)) * Re{ int e^{-i u x} phi(u + i a) / (i u - a) du }
% with R_a = 1 if a < 0 (left shift), R_a = 0 if a > 0 (right shift).
    int_kernel = @(csi) phi_cond(csi + 1i*a) ./ (1i*csi - a);
    fk_raw = arrayfun(int_kernel, g.xk);
    fk_raw(~isfinite(fk_raw)) = 0;
    fk = fk_raw .* exp(-1i * g.z1 * g.dx .* g.j);
 
    f_hat = g.dx .* exp(-1i * g.x1 * g.zk) .* fft(fk);
    f_hat = interp1(g.zk, f_hat, x, 'spline');
 
    P = Ra - real(f_hat(:) / (2*pi)) .* exp(a * x);
end