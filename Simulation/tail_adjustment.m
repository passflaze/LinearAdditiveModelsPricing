function [cdf_fine, x_fine] = tail_adjustment(x_grid, cdf_raw, refinement_factor)
%TAIL_ADJUSTMENT Cleans, refines, and extrapolates a numerically generated CDF.
%   [CDF_FINE, X_FINE] = TAIL_ADJUSTMENT(X_GRID, CDF_RAW, REFINEMENT_FACTOR)
%   takes a raw Cumulative Distribution Function (typically generated via FFT, 
%   which may contain numerical noise or oscillations in the tails) and 
%   reconstructs a clean, strictly monotonic CDF on a refined spatial grid.
%
%   The algorithm performs the following steps:
%   1. Scans outward from the midpoint to isolate the stable, monotonic core.
%   2. Interpolates the core onto a finer grid (defaults to 'spline', but falls 
%      back to 'pchip' to strictly preserve monotonicity if oscillations occur).
%   3. Estimates exponential decay rates (lambda) at the boundaries of the core.
%   4. Splices exponentially decaying tails to replace the numerical noise.
%   5. Enforces exact probability bounds [0, 1] at the extreme grid points.
%
%   Inputs:
%       x_grid            - Original spatial grid vector.
%       cdf_raw           - Raw CDF values evaluated on x_grid.
%       refinement_factor - (Optional) Integer multiplier to increase the grid 
%                           density. Default is 10.
%
%   Outputs:
%       cdf_fine          - The cleaned, strictly monotonic CDF.
%       x_fine            - The newly refined spatial grid.

    if nargin < 3
        refinement_factor = 10;
    end
    
    x_grid = x_grid(:);
    cdf_raw = cdf_raw(:);
    N = length(cdf_raw);
    midpoint = floor(N / 2);
    
    % =========================================================================
    % 2 & 3. OUTWARD SCAN 
    % =========================================================================
    % Isolate the strictly monotonic inner core of the CDF
    idx_b = midpoint;
    while idx_b > 2 && cdf_raw(idx_b - 1) > 0 && cdf_raw(idx_b - 1) < cdf_raw(idx_b)
        idx_b = idx_b - 1;
    end
    
    idx_e = midpoint;
    while idx_e < N - 1 && cdf_raw(idx_e + 1) < 1 && cdf_raw(idx_e + 1) > cdf_raw(idx_e)
        idx_e = idx_e + 1;
    end
    
    x_core = x_grid(idx_b:idx_e);
    cdf_core = cdf_raw(idx_b:idx_e);
    
    % =========================================================================
    % 4 & 5. GRID REFINEMENT & INTERPOLATION
    % =========================================================================
    step_raw = x_grid(2) - x_grid(1);
    step_fine = step_raw / refinement_factor;
    x_fine = (x_grid(1):step_fine:x_grid(end))';
    cdf_fine = zeros(size(x_fine));
    
    core_mask = (x_fine >= x_core(1)) & (x_fine <= x_core(end));
    
    cdf_test_spline = interp1(x_core, cdf_core, x_fine(core_mask), 'spline');
    
    inversions_spline = sum(diff(cdf_test_spline) < 0);
    if inversions_spline > 0
        fprintf('2. [WARNING] Spline interpolation caused %d monotonicity violations (oscillations).\n', inversions_spline);
        fprintf('   -> Automatic fallback to ''pchip'' applied to enforce monotonicity.\n');
        cdf_fine(core_mask) = interp1(x_core, cdf_core, x_fine(core_mask), 'pchip');
    else
        fprintf('2. Spline interpolation successful: No monotonicity violations detected.\n');
        cdf_fine(core_mask) = cdf_test_spline;
    end
    
    % =========================================================================
    % 6. LAMBDA CALCULATION 
    % =========================================================================
    xb = x_core(1); xb_next = x_core(2);
    Pb = cdf_core(1); Pb_next = cdf_core(2);
    
    xe = x_core(end); xe_prev = x_core(end-1);
    Pe = cdf_core(end); Pe_prev = cdf_core(end-1);
    
    lambda_minus = (log(Pb_next) - log(Pb)) / (xb_next - xb);
    lambda_plus  = (log(1 - Pe_prev) - log(1 - Pe)) / (xe - xe_prev);
    
    % =========================================================================
    % 7. EXPONENTIAL SPLICING & EXACT BOUNDARY FORCING
    % =========================================================================
    % Left Tail
    left_mask = (x_fine < xb);
    tail_sx_vals = Pb .* exp(lambda_minus .* (x_fine(left_mask) - xb));
    cdf_fine(left_mask) = tail_sx_vals;
    
    % Right Tail
    right_mask = (x_fine > xe);
    tail_dx_vals = 1 - (1 - Pe) .* exp(-lambda_plus .* (x_fine(right_mask) - xe));
    cdf_fine(right_mask) = tail_dx_vals;
    
    % Force bounds to exactly 0 and 1
    cdf_fine(1) = 0;
    cdf_fine(end) = 1;
    
end