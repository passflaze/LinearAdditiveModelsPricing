function [cdf_ext, x_ext] = smart_cdf_extrapolation(x_grid, cdf_raw, eps_tail, refinement_factor)
%SMART_CDF_EXTRAPOLATION  Extends a numerically-inverted CDF BEYOND the finite
%   FFT range using the analytic exponential tails (cf. point 0 and [5]).
%
%   The Lewis-FFT inversion returns the CDF only on the finite spatial grid
%   [x_grid(1), x_grid(end)]. Outside that range the raw CDF is truncated
%   (clamped to 0/1), so the inverse-CDF sampler can never produce increments
%   beyond the grid: the tail mass is lost and option prices are biased low.
%
%   The three models of Project 2 (MA, GL, AB) all have EXPONENTIAL tails:
%       1 - CDF(x) ~ exp(-p_plus  * x)   as x -> +inf
%           CDF(x) ~ exp( p_minus * x)   as x -> -inf .
%   We estimate the local decay rates (p_minus, p_plus) at the edges of the
%   stable monotonic core and splice exponential tails that EXTEND the support
%   until the CDF reaches eps_tail / 1 - eps_tail. This is the "smart
%   extrapolation" of the CDF outside the FFT range: it lets a not-so-large
%   FFT grid recover the tail mass that plain truncation throws away.
%
%   Difference vs tail_adjustment: tail_adjustment only cleans/refines the
%   tails WITHIN the original grid range (same support); here the support is
%   genuinely ENLARGED with the analytic exponential decay.
%
% INPUTS:
%   x_grid            - spatial grid of the raw CDF (vector)
%   cdf_raw           - raw CDF values on x_grid (vector)
%   eps_tail          - (optional, default 1e-6) residual tail probability at
%                       the extended edges
%   refinement_factor - (optional, default 5) grid-density multiplier
%
% OUTPUTS:
%   cdf_ext           - extended, strictly monotonic CDF (0 -> 1)
%   x_ext             - extended spatial grid (wider than x_grid)

    if nargin < 3 || isempty(eps_tail),          eps_tail          = 1e-6; end
    if nargin < 4 || isempty(refinement_factor), refinement_factor = 5;    end

    x_grid  = x_grid(:);
    cdf_raw = min(max(cdf_raw(:), 0), 1);
    N       = numel(cdf_raw);
    mid     = floor(N / 2);

    % --- Isolate the strictly-monotonic, strictly-interior core --------------
    idx_b = mid;
    while idx_b > 2 && cdf_raw(idx_b-1) > 0 && cdf_raw(idx_b-1) < cdf_raw(idx_b)
        idx_b = idx_b - 1;
    end
    idx_e = mid;
    while idx_e < N-1 && cdf_raw(idx_e+1) < 1 && cdf_raw(idx_e+1) > cdf_raw(idx_e)
        idx_e = idx_e + 1;
    end

    x_core   = x_grid(idx_b:idx_e);
    cdf_core = cdf_raw(idx_b:idx_e);

    xb = x_core(1);    Pb = cdf_core(1);
    xe = x_core(end);  Pe = cdf_core(end);

    % --- Local exponential decay rates at the core edges ---------------------
    lambda_minus = (log(cdf_core(2))   - log(Pb))   / (x_core(2)   - xb);   % left  tail
    lambda_plus  = (log(1-cdf_core(end-1)) - log(1-Pe)) / (xe - x_core(end-1)); % right tail

    % --- Extension lengths so the spliced tails reach eps_tail ---------------
    ext_left  = 0;
    if isfinite(lambda_minus) && lambda_minus > 0
        ext_left = max(0, (log(Pb) - log(eps_tail)) / lambda_minus);
    end
    ext_right = 0;
    if isfinite(lambda_plus) && lambda_plus > 0
        ext_right = max(0, (log(1-Pe) - log(eps_tail)) / lambda_plus);
    end

    % --- Refined, ENLARGED grid ---------------------------------------------
    step_fine = (x_grid(2) - x_grid(1)) / refinement_factor;
    x_ext     = ((xb - ext_left) : step_fine : (xe + ext_right))';
    cdf_ext   = zeros(size(x_ext));

    core_mask  = (x_ext >= xb) & (x_ext <= xe);
    left_mask  = (x_ext <  xb);
    right_mask = (x_ext >  xe);

    % PCHIP keeps the core monotone (no spline overshoot)
    cdf_ext(core_mask)  = interp1(x_core, cdf_core, x_ext(core_mask), 'pchip');
    cdf_ext(left_mask)  = Pb .* exp(lambda_minus .* (x_ext(left_mask) - xb));
    cdf_ext(right_mask) = 1 - (1 - Pe) .* exp(-lambda_plus .* (x_ext(right_mask) - xe));

    % --- Enforce monotonicity and exact [0,1] bounds -------------------------
    cdf_ext      = min(max(cdf_ext, 0), 1);
    cdf_ext      = cummax(cdf_ext);
    cdf_ext(1)   = 0;
    cdf_ext(end) = 1;
end
