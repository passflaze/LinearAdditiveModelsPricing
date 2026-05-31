function call_GL_final = price_GL(alpha, beta, M, dz, discount_factors, sigma_ATM, yf, modified_moneyness)
% PRICE_GL European call prices under the Generalized Logistic (GL) model
% via a Lewis-Bachelier FFT with **double damping**.
%
% A single damping shift (e.g. shift = beta/2, lower contour) is appropriate
% for the OTM call region but produces an exp(-beta*k/2) prefactor that
% diverges for k << 0 (ITM call / OTM put region). To keep numerical
% accuracy across the full surface we run two Lewis FFTs:
%
%   * Lower contour, shift_call = +beta/2, used for k_tilde >= 0:
%       C(k_tilde) = -(1/(2*pi)) * exp(-beta*k_tilde/2) *
%                    integral over Im(w)=-beta/2 of e^{-i w k_tilde} phi(w)/w^2 dw
%
%   * Upper contour, shift_put = -alpha/2, used for k_tilde < 0:
%       P(k_tilde) = -(1/(2*pi)) * exp(+alpha*k_tilde/2) *
%                    integral over Im(w)=+alpha/2 of e^{-i w k_tilde} phi(w)/w^2 dw
%
% The call price for ITM strikes is then recovered via Bachelier put-call
% parity in zeta-space (martingale, E[zeta]=0):
%       C_zeta(k_tilde) = -k_tilde + P_zeta(k_tilde),
% which in market units becomes
%       C(K) = B * (F - K) + P(K).
%
% INPUTS:
%   alpha, beta        : GL shape parameters (>0). Analyticity strip of phi
%                        is Im(z) in (-beta, alpha), so the two shifts
%                        +beta/2 and -alpha/2 both lie strictly inside.
%   M                  : N_FFT = 2^M
%   dz                 : output (k_tilde) grid step
%   discount_factors   : (T x 1) discount factors
%   sigma_ATM          : (T x 1) Bachelier ATM vols
%   yf                 : (T x 1) year fractions
%   modified_moneyness : (T x K) market chi = (K - F)/(sigma_ATM*sqrt(t))
%
% OUTPUT:
%   call_GL_final      : (T x K) model call prices on the same grid as the
%                        input moneyness (NaN where moneyness is NaN).

    % =========================================================================
    % STEP 0: NORMALIZATION CONSTANT
    % =========================================================================
    I0 = I0_GL(alpha, beta);
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
