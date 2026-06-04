function call_GL_final = price_GL(alpha, beta, M, dz, discount_factors, sigma_ATM, yf, modified_moneyness)
% PRICE_GL  European call surface under the Generalized Logistic model via a
%   double-damping Lewis-FFT. Two contours keep accuracy across the surface:
%     k_tilde >= 0 (OTM): lower contour shift +beta/2;
%     k_tilde <  0 (ITM): upper contour shift -alpha/2, plus parity
%                         C = B*(F-K) + P.
%   phi strip is Im(z) in (-beta, alpha), so both shifts lie strictly inside.
%
% INPUTS:
%   alpha, beta        : GL shape parameters (> 0)
%   M, dz              : FFT exponent (N = 2^M) and output grid step
%   discount_factors, sigma_ATM, yf : (T x 1) per maturity
%   modified_moneyness : (T x K) chi = (K-F)/(sigma_ATM*sqrt(t)), NaN if no quote
% OUTPUT:
%   call_GL_final      : (T x K) call prices (NaN where chi is NaN)

    % =========================================================================
    % STEP 0: NORMALIZATION CONSTANT
    % =========================================================================
    I0 = I0_GL([alpha; beta]);
    if isnan(I0) || isinf(I0) || I0 == 0
        error('PricingEngine:I0_Invalid', ...
              'I0 evaluated to an unphysical value (I0 = %f).', I0);
    end

    shift_call = beta  / 2;   % positive shift -> lower contour
    shift_put  = -alpha / 2;  % negative shift -> upper contour

    % =========================================================================
    % STEP 1: FFT GRID (shared between the two transforms)
    % =========================================================================
    N  = 2^M;
    dx = (2*pi) / (N * dz);

    zn = (dz * (N-1)) / 2;
    z1 = -zn;
    z_grid = z1 : dz : zn;

    xn = (dx * (N-1)) / 2;
    x1 = -xn;
    x_grid = x1 : dx : xn;

    j_arr     = 0:N-1;
    prefactor = dx * exp(-1i * x1 * z_grid);
    twist     = exp(-1i * z1 * dx * j_arr);

    % =========================================================================
    % STEP 2: TWO LEWIS FFTs (one per damping)
    % =========================================================================
    raw_call = lewis_fft_(x_grid, alpha, beta, shift_call, prefactor, twist);
    raw_put  = lewis_fft_(x_grid, alpha, beta, shift_put,  prefactor, twist);

    src_call = real(raw_call);
    src_put  = real(raw_put);

    % =========================================================================
    % STEP 3: BRANCH SELECTION ON k_tilde = I0 * chi
    % =========================================================================
    k_tilde   = modified_moneyness * I0;
    mask_data = ~isnan(modified_moneyness);
    mask_call = mask_data & (k_tilde >= 0);
    mask_put  = mask_data & (k_tilde <  0);

    interp_call = nan(size(modified_moneyness));
    interp_put  = nan(size(modified_moneyness));

    ok_call = isfinite(src_call);
    ok_put  = isfinite(src_put);

    if any(mask_call(:)) && nnz(ok_call) >= 2
        interp_call(mask_call) = interp1(z_grid(ok_call), src_call(ok_call), ...
                                         k_tilde(mask_call), 'spline');
    end
    if any(mask_put(:))  && nnz(ok_put)  >= 2
        interp_put(mask_put)   = interp1(z_grid(ok_put),  src_put(ok_put), ...
                                         k_tilde(mask_put),  'spline');
    end

    % =========================================================================
    % STEP 4: ASSEMBLE THE CALL PRICE
    % =========================================================================
    % Maturity-dependent factor: B * sigma_GL * sqrt(t) with sigma_GL = sigma_ATM/I0.
    B_sigma_sqrt = discount_factors .* (sigma_ATM / I0) .* sqrt(yf);   % (T x 1)

    pref_call = -B_sigma_sqrt .* exp(-shift_call .* k_tilde) / (2*pi);
    pref_put  = -B_sigma_sqrt .* exp(-shift_put  .* k_tilde) / (2*pi);
    intrinsic = -B_sigma_sqrt .* k_tilde;   % equals B * (F - K) by construction

    call_GL_final = nan(size(modified_moneyness));
    call_GL_final(mask_call) = pref_call(mask_call) .* interp_call(mask_call);
    call_GL_final(mask_put)  = intrinsic(mask_put) + ...
                               pref_put(mask_put) .* interp_put(mask_put);

end

% -------------------------------------------------------------------------
function out = lewis_fft_(x_grid, alpha, beta, shift, prefactor, twist)
% Evaluate dx * sum_j integrand(x_j) * exp(-i x_j z_l) on the z_grid implied
% by prefactor/twist. integrand(x) = phi(x - i*shift)/(x - i*shift)^2.
    cf_vals   = cf_GL_cal(x_grid - 1i*shift, alpha, beta);
    integrand = cf_vals ./ ((x_grid - 1i*shift).^2);
    out       = prefactor .* fft(integrand .* twist);
end
