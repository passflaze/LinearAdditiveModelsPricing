function call_GL_final = price_GL(alpha, beta, M, dz, discount_factors, sigma_ATM, yf, modified_moneyness)
% PRICE_GL Computes European Call option prices under the Generalized Laplace model via FFT.
%
% This function implements a Fast Fourier Transform (FFT) approach (similar to 
% Carr-Madan / Lewis) to price options by evaluating the characteristic function 
% on a discretized frequency grid and interpolating the results back to the market grid.
%
% INPUTS:
%   alpha              : Left shape parameter of the GL distribution
%   beta               : Right shape parameter of the GL distribution (used for the damping shift)
%   M                  : Power of 2 determining the number of FFT grid points (N = 2^M)
%   dz                 : Step size for the spatial log-moneyness grid (z_grid)
%   discount_factors   : Discount factor for the given maturity
%   sigma_ATM          : At-The-Money implied volatility parameter
%   yf                 : Year fraction (time to maturity)
%   modified_moneyness : Matrix or vector of modified moneyness used in the prefactor scaling
%   x_target           : Target market modified moneyness points where the final prices are interpolated
%
% OUTPUT:
%   call_GL_interp     : Interpolated Call option prices at the requested x_target points

    % =========================================================================
    % STEP 1: CALCULATE I0 (Normalization Constant) 
    % =========================================================================
    shift = beta / 2;             
    
    % =========================================================================
    % STEP 1: CALCULATE I0 (Normalization Constant)
    % =========================================================================
    % Convention from Baviera & Massaria (2026), Eq. 14-15: I0 := sqrt(2*pi) * E[zeta_+]
    % so that sigma_t = sigma_ATM / I0 enforces the ATM Bachelier condition exactly.
    % Without the sqrt(2*pi) factor the whole price surface is off by ~2.5x and the
    % interpolation point modified_moneyness*I0 lands at the wrong moneyness.
    integrand_mean = @(x) pdf_GL(alpha, beta, x) .* x;
    I0 = sqrt(2*pi) * quadgk(integrand_mean, 0, inf);

    % Check for mathematical invalidity of the constant
    if isnan(I0) || isinf(I0) || I0 == 0
        error('PricingEngine:I0_Invalid', 'I0 evaluated to an unphysical value (I0 = %f). Check cf_GL parameters.', I0);
    end
    
    % =========================================================================
    % STEP 2: FFT GRID SETUP 
    % =========================================================================
    N = 2^(M);  
    dx = (2*pi) / (N * dz);         
    
    zn = (dz * (N-1)) / 2;
    z1 = -zn;
    z_grid = z1 : dz : zn;
    
    xn = (dx * (N-1)) / 2;
    x1 = -xn;
    x_grid = x1 : dx : xn;
    
    % =========================================================================
    % STEP 3: CALCULATE CALL PRICES VIA FFT
    % =========================================================================
    prefactor = dx * exp(-1i * x1 * (z_grid));
    preprefactor = discount_factors .* (sigma_ATM / I0) .* sqrt(yf) .* (-exp(-shift .* modified_moneyness * I0) / (2*pi));
    
    fourier_function1 = cf_GL(x_grid - 1i*shift, alpha, beta) ;
    fourier_function = fourier_function1./ ((x_grid - 1i*shift).^2);
    
    j_minus_1 = 0:N-1;
    input_fft = fourier_function .* exp(-1i * z1 * dx * j_minus_1);
    
    fft_call_prices = fft(input_fft);
    call_GL = prefactor .* fft_call_prices;
    
    % =========================================================================
    % STEP 4: INTERPOLATION TO MARKET GRID
    % =========================================================================
    src = real(call_GL);
    src_ok = isfinite(src);                              % filtra sorgente
    mask = ~isnan(modified_moneyness);
    call_GL_interp = nan(size(modified_moneyness));
    
    if nnz(src_ok) >= 2 && any(mask(:))
        call_GL_interp(mask) = interp1(z_grid(src_ok), src(src_ok), ...
                                        modified_moneyness(mask)*I0, 'spline');
    end
    
    call_GL_final = preprefactor .* call_GL_interp;
end

