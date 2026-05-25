function Z = sample_from_cdf(x_grid, cdf, Nsim)
% Sample Nsim increments via inverse-CDF (probability integral transform).
%
% Implements Baviera-Manzoni (2026) Algorithm 1, steps 2-4:
%   2. Restrict x-grid to the largest monotone-in-[0,1] segment.
%   3. Cubic spline interpolation (not-a-knot, MATLAB default).
%      Exponential tail extrapolation outside [xb, xe].
%   4. Draw U ~ Uniform(0,1) and invert: Z = F^{-1}(U).
%
% Inputs:
%   x_grid : N x 1 vector of moneyness points (output of ccdf_AB_FFT)
%   cdf    : N x 1 CDF values at x_grid
%   Nsim   : number of samples
%
% Output:
%   Z      : Nsim x 1 vector of simulated increments

x_grid = x_grid(:);
cdf    = cdf(:);

% restrict to monotone segment in [0, 1] 
in_bounds = (cdf >= 0) & (cdf <= 1);
idx_b = find(in_bounds, 1, 'first');
idx_e = find(in_bounds, 1, 'last');

x_v = x_grid(idx_b:idx_e);
c_v = cdf(idx_b:idx_e);

% enforce strict monotonicity (remove plateaus and micro-dips)
c_v  = cummax(c_v);
keep = [true; diff(c_v) > 0];
x_v  = x_v(keep);
c_v  = c_v(keep);

% rule-of-thumb check (Baviera-Manzoni 2026, Sec. 5.3)
tail_err = max(c_v(1), 1 - c_v(end));
if tail_err > 1e-4
    warning('sample_from_cdf:tails', ...
        'Tail mass = %.2e > 1e-4. Widen x_grid or refine CDF grid.', tail_err);
end

% Step 3 + 4: sample via inverse spline 
U = rand(Nsim, 1);
Z = interp1(c_v, x_v, U, 'spline');

% Exponential tail extrapolation for U outside [c_v(1), c_v(end)] --
% Left tail:  F(x) = c_v(1) * exp(lam_m * (x - xb))  for x <= xb
% Right tail: 1-F(x) = (1-c_v(end)) * exp(-lam_p * (x - xe))  for x >= xe
xb = x_v(1);   xe = x_v(end);
Pb = c_v(1);   Pe = c_v(end);

lam_m = (log(c_v(2))     - log(Pb))       / (x_v(2)     - xb);
lam_p = (log(1 - c_v(end-1)) - log(1 - Pe)) / (xe - x_v(end-1));  % > 0

left  = U < Pb;
right = U > Pe;
Z(left)  = xb + log(U(left)  / Pb)       / lam_m;
Z(right) = xe - log((1 - U(right)) / (1 - Pe)) / lam_p;

Z = Z(:);
end
