function cdf = cdf_AB_increments(eta, k, T1, T2, sigma_T1, sigma_T2, x, fwd_factor)
    % Conditional CDF F_{T2|T1}(x) of the AB log-price increment between
    % reset dates T1 and T2, via Lewis-FFT.
    %
    % Two-shift reconstruction (Baviera-Manzoni 2026, Sec. 5.3): computes the
    % CDF twice, once with a shift on the negative side of the analyticity
    % strip (accurate on the LEFT tail) and once on the positive side
    % (accurate on the RIGHT tail), then glues the two at the median (x = 0).
    %
    % sigma_T1, sigma_T2 = entries of sigma_t from calibrateAB (sigmaATM/I_0).
    %   [~, ~, sigma_t, ~] = calibrateAB(...);
    %   cdf = cdf_AB_increments(eta, k, T1, T2, sigma_t(iT1), sigma_t(iT2), x);
    %
    % x: column vector of log-price increment points where to evaluate the CDF.
    %
    % fwd_factor : (optional, default 1) Lemma 2 rescaling (Forward.pdf).
    %   Under Lemma 2 strict, f_{T1,T2} = fwd_factor * f_{T1,T1} with
    %   fwd_factor = B(0,T1)/B(0,T2). The increment W = f_{T2,T2} - f_{T1,T2}
    %   of the AB-T2 process has CF phi_T2(u) / phi_T1(fwd_factor * u).

    if nargin < 8 || isempty(fwd_factor)
        fwd_factor = 1;
    end

    phi_T2   = charateristic_function_AB(T2, k, eta, sigma_T2);
    if T1 == 0
        phi_cond = phi_T2;
    else
        phi_T1   = charateristic_function_AB(T1, k, eta, sigma_T1);
        phi_cond = @(u) phi_T2(u) ./ phi_T1(fwd_factor * u);
    end

    % FFT grid (single wide "Lewis recipe", flag kept only for signature compat).
    % For fixed N the spatial range and Fourier step are LOCKED:
    %   spatial range = N*dz ,   dx = 2*pi/(N*dz) = 2*pi/range .
    % The conditional increment W is in dollars (sigma_W ~ 10$), so an
    % evaluation grid at +/-10 sigma reaches ~+/-106$. The old flag=0 grid
    % (dz=0.0025) gave a spatial range of only +/-82$: interp1 then EXTRAPOLATED
    % the right tail (biasing E[(W)+] ~10% low vs the Lewis reference) and its
    % dx=0.038 under-resolved the integrand near-pole at u=0 (width ~|a|~0.03).
    % dz=0.05 -> range +/-1638$ (covers any dollar grid) and dx~0.0019 (resolves
    % the peak); identical to the validated lewis_FFT_AB grid.
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

    % Analyticity strip of the conditional CF
    % phi_T2 strip in u: |Im u| < p+/(sigma_T2*sqrt(T2)).
    % phi_T1(fwd_factor*u) strip in u: |Im u| < p+/(fwd_factor*sigma_T1*sqrt(T1)).
    % Under Lemma 2, take the binding (tightest) side.
    p_plus  = eta + sqrt(eta^2 + 1/k);
    p_minus = -eta + sqrt(eta^2 + 1/k);
    sT2     = sigma_T2 * sqrt(T2);

    if T1 == 0
        strip_plus  = p_plus  / sT2;
        strip_minus = p_minus / sT2;
    else
        sT1_eff     = fwd_factor * sigma_T1 * sqrt(T1);
        strip_plus  = min(p_plus  / sT2, p_plus  / sT1_eff);
        strip_minus = min(p_minus / sT2, p_minus / sT1_eff);
    end

    a_neg = -0.49 * strip_plus;    % left edge of strip:  good for RIGHT tail, Ra = 1
    a_pos = +0.49 * strip_minus;   % right edge of strip: good for LEFT tail,  Ra = 0

    % Two FFT reconstructions
    x = x(:);
    cdf_right = one_shift(phi_cond, x, a_neg, 1, grid);   % accurate for x > 0
    cdf_left  = one_shift(phi_cond, x, a_pos, 0, grid);   % accurate for x < 0

    % Glue with a SMOOTH blend, not a hard switch at x = 0.
    % A brusque switch leaves a downward step at the median: the raw CDF is no
    % longer monotone, so trapz(1-cdf) (FFT reference) and the cummax-cleaned
    % CDF used by sample_from_cdf integrate different functions -> the MC
    % forward-start price drifts ~10% above the FFT reference.
    % Both reconstructions are accurate everywhere except their own far tail,
    % so a wide tanh blend (scale ~ a fraction of the x-range) suppresses each
    % side's oscillation by averaging where both are good, and recovers the
    % accurate tail far from 0 (w->1 right, w->0 left). The result is C^1 and
    % monotone, so cummax in sample_from_cdf becomes a no-op.
    scale = (max(x) - min(x)) / 20;          % gentle transition (~1 std)
    w     = 0.5 * (1 + tanh(x / scale));     % 0 for x<<0, 1 for x>>0
    cdf   = (1 - w) .* cdf_left + w .* cdf_right;

    % numerical hygiene: clip to [0,1] and enforce monotonicity
    cdf = min(max(cdf, 0), 1);
    cdf = cummax(cdf);

end

% ----------------------------------------------------------------------
function P = one_shift(phi_cond, x, a, Ra, g)
% Single FFT reconstruction of the CDF using Baviera-Manzoni eq. (13)-(15):
%   P(x) = Ra - (e^{-ap x}/(2 pi)) * int e^{-iux} phi(u - i ap) / (i (u - i ap)) du
% In code convention ( a = -ap ):
%   P(x) = Ra - real(FFT_output / (2 pi)) * exp(a * x)
% with Ra = 1 for a < 0 and Ra = 0 for a > 0.
    % vectorised CF call: phi_cond returns column N x 1; reshape to row for FFT
    phi_vals = phi_cond(g.xk + 1i*a);
    phi_vals = phi_vals(:).';
    fk_raw = phi_vals ./ (1i*g.xk - a);
    fk_raw(~isfinite(fk_raw)) = 0;
    fk = fk_raw .* exp(-1i * g.z1 * g.dx .* g.j);

    f_hat = g.dx .* exp(-1i * g.x1 * g.zk) .* fft(fk);
    f_hat = interp1(g.zk, f_hat, x, 'spline');

    P = Ra - real(f_hat(:) / (2*pi)) .* exp(a * x);
end
