function cdf = ccdf_increment_FFT(model, params, T1, T2, sigma_T1, sigma_T2, x, fwd_factor)
% CCDF_INCREMENT_FFT  Conditional CDF F_{T2|T1}(x) of the Linear Additive
% log-price increment W = f_{T2,T2} - f_{T1,T2} between reset dates T1 -> T2,
% via Lewis-FFT. Model-agnostic: works for 'AB' and 'GL' (the only model
% dependence is the marginal CF / strip returned by model_marginal_cf).
%
% Model-agnostic two-shift CDF reconstruction:
%   * two-shift reconstruction (Baviera-Manzoni 2026, Sec.5.3): the CDF is
%     built twice -- once shifting the contour on the NEGATIVE side of the
%     analyticity strip (accurate on the RIGHT tail) and once on the POSITIVE
%     side (accurate on the LEFT tail);
%   * the two halves are spliced with a SMOOTH tanh blend (not a hard switch
%     at x=0) so the result is C^1 and monotone -> no spurious step at the
%     median that would bias the survival-integral price;
%   * Lemma 2 forward rescaling (Forward.pdf): under f_{T1,T2} = fwd_factor *
%     f_{T1,T1} with fwd_factor = B(0,T1)/B(0,T2), the conditional CF is
%       phi_W(u) = phi_T2(u) / phi_T1(fwd_factor * u).
%
% INPUTS
%   model      : 'AB' or 'GL'
%   params     : model parameter column vector (see model_marginal_cf)
%   T1, T2     : reset and maturity (year fractions). T1 = 0 -> marginal CDF.
%   sigma_T1   : scale at T1 (pass 0 if T1 == 0).
%   sigma_T2   : scale at T2.
%   x          : column vector of $-increments where to evaluate the CDF.
%   fwd_factor : (optional, default 1) Lemma 2 rescaling B(0,T1)/B(0,T2).
%
% OUTPUT
%   cdf        : column vector, F_W(x), clipped to [0,1] and monotone.

    if nargin < 8 || isempty(fwd_factor)
        fwd_factor = 1;
    end

    % --- marginal CFs + analyticity-strip constants -----------------------
    [phi_T2, c_neg, c_pos] = model_marginal_cf(model, params, T2, sigma_T2);

    if T1 == 0
        phi_cond = phi_T2;                              % marginal at T2
    else
        phi_T1   = model_marginal_cf(model, params, T1, sigma_T1);
        phi_cond = @(u) phi_T2(u) ./ phi_T1(fwd_factor * u);
    end

    % --- FFT grid (single wide "Lewis recipe", dollar increments) ---------
    % For fixed N the spatial range and Fourier step are LOCKED:
    %   spatial range = N*dz ,  dx = 2*pi/(N*dz) .
    % dz=0.05 -> range ~+/-1638$ (covers any dollar grid) and dx~0.0019
    % (resolves the near-pole at u=0); identical to the validated AB/GL grids.
    M  = 16;
    dz = 0.05;
    N  = 2^M;
    dx = 2*pi / (N*dz);
    z1 = -dz * (N-1) / 2;
    x1 = -dx * (N-1) / 2;
    j  = 0:N-1;
    zk = z1 + dz*j;
    xk = x1 + dx*j;
    grid = struct('zk', zk, 'xk', xk, 'x1', x1, 'z1', z1, 'dx', dx, 'j', j);

    % --- binding analyticity strip ----------------------------------------
    % phi_T2 strip in u: |Im u| < c/(sigma_T2*sqrt(T2)).
    % phi_T1(fwd*u) strip in u: |Im u| < c/(fwd*sigma_T1*sqrt(T1)).
    % Under Lemma 2 take the tightest (min) side.
    sT2 = sigma_T2 * sqrt(T2);
    if T1 == 0
        strip_neg = c_neg / sT2;
        strip_pos = c_pos / sT2;
    else
        sT1_eff   = fwd_factor * sigma_T1 * sqrt(T1);
        strip_neg = min(c_neg / sT2, c_neg / sT1_eff);
        strip_pos = min(c_pos / sT2, c_pos / sT1_eff);
    end

    a_neg = -0.49 * strip_neg;    % left edge of strip:  good for RIGHT tail, Ra = 1
    a_pos = +0.49 * strip_pos;    % right edge of strip: good for LEFT tail,  Ra = 0

    % --- two FFT reconstructions ------------------------------------------
    x = x(:);
    cdf_right = one_shift(phi_cond, x, a_neg, 1, grid);   % accurate for x > 0
    cdf_left  = one_shift(phi_cond, x, a_pos, 0, grid);   % accurate for x < 0

    % --- smooth tanh blend (C^1, monotone) --------------------------------
    scale = (max(x) - min(x)) / 20;          % gentle transition (~1 std)
    w     = 0.5 * (1 + tanh(x / scale));     % 0 for x<<0, 1 for x>>0
    cdf   = (1 - w) .* cdf_left + w .* cdf_right;

    % numerical hygiene: clip to [0,1] and enforce monotonicity
    cdf = min(max(cdf, 0), 1);
    cdf = cummax(cdf);
end

% ----------------------------------------------------------------------
function P = one_shift(phi_cond, x, a, Ra, g)
% Single FFT reconstruction of the CDF (Baviera-Manzoni eq.13-15):
%   P(x) = Ra - real(FFT_output / (2 pi)) * exp(a * x)
% with Ra = 1 for a < 0 (left shift) and Ra = 0 for a > 0 (right shift).
    phi_vals = phi_cond(g.xk + 1i*a);
    phi_vals = phi_vals(:).';                 % force row for FFT
    fk_raw   = phi_vals ./ (1i*g.xk - a);
    fk_raw(~isfinite(fk_raw)) = 0;
    fk = fk_raw .* exp(-1i * g.z1 * g.dx .* g.j);

    f_hat = g.dx .* exp(-1i * g.x1 * g.zk) .* fft(fk);
    f_hat = interp1(g.zk, f_hat, x, 'spline');

    P = Ra - real(f_hat(:) / (2*pi)) .* exp(a * x);
end
