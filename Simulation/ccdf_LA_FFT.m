function cdf = ccdf_LA_FFT(spec, T1, T2, sigma_T1, sigma_T2, x, fwd_factor)
%CCDF_LA_FFT  Conditional CDF F_{T2|T1}(x) of a Linear Additive log-price
%   increment between reset dates T1 and T2, via the Lewis-FFT two-shift
%   reconstruction. Unified engine for the AB and GL models: it is a direct
%   generalization of the validated ccdf_AB_FFT (same Lewis recipe, same smooth
%   tanh gluing, same Lemma-2 forward rescaling) where the only model-specific
%   inputs -- the marginal CF and the analyticity-strip tail rates -- are read
%   from the SPEC struct (see la_model_spec).
%
%   Two-shift reconstruction (Baviera-Manzoni 2026, Sec. 5.3): the CDF is
%   computed twice, once with the contour shifted toward the lower edge of the
%   analyticity strip (accurate on the RIGHT tail) and once toward the upper edge
%   (accurate on the LEFT tail), then blended with a wide tanh window so the
%   result is C^1, monotone and free of the median step that a hard switch leaves.
%
%   INPUTS
%     spec       : model spec struct from la_model_spec (fields .cf, .lamL, .lamR)
%     T1, T2     : reset and maturity year fractions (T1 = 0 -> marginal CDF of f_T2)
%     sigma_T1   : model scale at T1 (= sigma_ATM(T1)/I0); pass 0 if T1 == 0
%     sigma_T2   : model scale at T2 (= sigma_ATM(T2)/I0)
%     x          : column vector of increment points where to evaluate the CDF
%     fwd_factor : (optional, default 1) Lemma-2 rescaling (Forward.pdf):
%                  f_{T1,T2} = fwd_factor * f_{T1,T1} with fwd_factor = B(0,T1)/B(0,T2).
%                  The increment W = f_{T2,T2} - f_{T1,T2} then has conditional CF
%                  phi_cond(u) = phi_T2(u) / phi_T1(fwd_factor * u).
%
%   OUTPUT
%     cdf : column vector of CDF values at x (clipped to [0,1], monotone).

    if nargin < 7 || isempty(fwd_factor)
        fwd_factor = 1;
    end

    % --- conditional characteristic function (additivity + Lemma 2) -----------
    phi_T2 = @(u) spec.cf(u, sigma_T2, T2);
    if T1 == 0
        phi_cond = phi_T2;
        s_den    = sigma_T2 * sqrt(T2);
    else
        phi_T1   = @(u) spec.cf(u, sigma_T1, T1);
        phi_cond = @(u) phi_T2(u) ./ phi_T1(fwd_factor * u);
        % Binding (tightest) analyticity strip across the two factors: phi_T2
        % constrains a via sigma_T2*sqrt(T2); phi_T1(fwd_factor*u) via
        % fwd_factor*sigma_T1*sqrt(T1). Both edges scale as 1/max(.) .
        s_den = max(sigma_T2 * sqrt(T2), fwd_factor * sigma_T1 * sqrt(T1));
    end

    % --- analyticity strip in the contour shift a ----------------------------
    %   density ~ exp(-lamR*z) (right) -> a > -lamR/s_den (lower edge)
    %   density ~ exp(+lamL*z) (left)  -> a < +lamL/s_den (upper edge)
    strip_lower = -spec.lamR / s_den;     % < 0
    strip_upper =  spec.lamL / s_den;     % > 0

    a_neg = 0.49 * strip_lower;   % lower-edge shift: accurate RIGHT tail, Ra = 1
    a_pos = 0.49 * strip_upper;   % upper-edge shift: accurate LEFT  tail, Ra = 0

    % --- FFT grid (Lewis recipe, locked range) -------------------------------
    % For fixed N the spatial range N*dz and Fourier step dx = 2*pi/(N*dz) are
    % locked. dz = 0.05 -> range +/-1638$ (covers any dollar grid, no interp1
    % extrapolation) and dx ~ 0.0019 (resolves the near-pole at u = 0). Identical
    % grid to the validated lewis_FFT_AB / ccdf_AB_FFT.
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

    % --- two reconstructions + smooth blend (as in ccdf_AB_FFT) --------------
    x = x(:);
    cdf_right = one_shift(phi_cond, x, a_neg, 1, grid);   % accurate for x > 0
    cdf_left  = one_shift(phi_cond, x, a_pos, 0, grid);   % accurate for x < 0

    scale = (max(x) - min(x)) / 20;          % gentle ~1-std transition
    w     = 0.5 * (1 + tanh(x / scale));     % 0 for x<<0, 1 for x>>0
    cdf   = (1 - w) .* cdf_left + w .* cdf_right;

    % numerical hygiene: clip to [0,1] and enforce monotonicity
    cdf = min(max(cdf, 0), 1);
    cdf = cummax(cdf);
end

% -------------------------------------------------------------------------
function P = one_shift(phi_cond, x, a, Ra, g)
% Single FFT reconstruction of the CDF (Baviera-Manzoni eq. 13-15):
%   P(x) = Ra - real( FFT_output / (2*pi) ) .* exp(a*x)
% with Ra = 1 for a < 0 (lower-edge shift) and Ra = 0 for a > 0 (upper-edge).
    phi_vals = phi_cond(g.xk + 1i*a);
    phi_vals = phi_vals(:).';                  % force row regardless of model
    fk_raw = phi_vals ./ (1i*g.xk - a);
    fk_raw(~isfinite(fk_raw)) = 0;
    fk = fk_raw .* exp(-1i * g.z1 * g.dx .* g.j);

    f_hat = g.dx .* exp(-1i * g.x1 * g.zk) .* fft(fk);
    f_hat = interp1(g.zk, f_hat, x, 'spline');

    P = Ra - real(f_hat(:) / (2*pi)) .* exp(a * x);
end
