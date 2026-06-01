function price = lewis_FFT_call(cf, M, dz, params, scalefactor, strikes, doubleshift, model)
% LEWIS_FFT_CALL_PROVA  Computes call prices via FFT inversion of the CF
%                       using the Lewis (2001) formula with dual damping shifts.
%
%
% INPUTS:
%   cf          - (function handle) CF of the increment: cf(u, params, scalefactor)
%   M           - (scalar) grid size exponent, N = 2^M (typical: 12-15)
%   dz          - (scalar) step size in the z-domain (log-strike grid)
%   params      - (vector) model parameters passed to cf
%   scalefactor - [scale_t1, scale_t2] scaling factors at t1 and t2
%   strikes     - (vector) log-strike values at which to interpolate prices
%   doubleshift - (logical) if true, uses both shifts and blends the grids
%   model       - (string) model identifier: 'MA', 'GL', or 'AB'
%
% OUTPUT:
%   price       - (vector) call prices interpolated at the requested strikes

    % FFT grid setup 
    N  = 2^M;
    dx = (2*pi) / (N * dz);

    zn     = (dz * (N-1)) / 2;
    z1     = -zn;
    z_grid = z1 : dz : zn;

    xn     = (dx * (N-1)) / 2;
    x1     = -xn;
    x_grid = x1 : dx : xn;

    % Damping shifts 
    switch model
        case 'MA'
            % params(1) = alpha, params(2) = beta
            shift_pos =  0.5* params(1) / scalefactor(2);
            shift_neg = -0.5 * params(2) / scalefactor(2);

        case 'GL'
            % params(1) = alpha, params(2) = beta
            shift_pos =  0.5 * params(1) / scalefactor(2);
            shift_neg = -0.5* params(2) / scalefactor(2);

        case 'AB'
            % params(1) = k (kappa), params(2) = eta
            kappa = params(1);
            eta   = params(2);
            shift_pos =  0.5 * (-eta + sqrt(eta^2 + 1/kappa)) / scalefactor(2);
            shift_neg = -0.5 * (eta + sqrt(eta^2 + 1/kappa)) / scalefactor(2);
    end

    % --- Negative shift: FFT inversion ---
    j_minus_1 = 0:N-1;
    prefactor        = dx * exp(-1i * x1 * z_grid);
    preprefactor_neg = -exp(shift_neg * z_grid) / (2*pi);

    x_grid_shifted_neg   = x_grid + 1i * shift_neg;
    fourier_function_neg = cf(x_grid_shifted_neg, params, scalefactor) ./ ...
                           (x_grid_shifted_neg.^2);
    input_fft_neg        = fourier_function_neg .* exp(-1i * z1 * dx * j_minus_1);

    price_clean_neg = preprefactor_neg .* real(prefactor .* fft(input_fft_neg));

    % --- Positive shift + blending (doubleshift mode) ---
    if doubleshift
        preprefactor_pos = -exp(shift_pos * z_grid) / (2*pi);

        x_grid_shifted_pos   = x_grid + 1i * shift_pos;
        fourier_function_pos = cf(x_grid_shifted_pos, params, scalefactor) ./ ...
                               (x_grid_shifted_pos.^2);
        input_fft_pos        = fourier_function_pos .* exp(-1i * z1 * dx * j_minus_1);

        price_clean_pos = -z_grid + preprefactor_pos .* real(prefactor .* fft(input_fft_pos));

        % Sigmoid blend in the transition region [-1, 1].
        % The positive shift (a>0, Ra=-z) is accurate on the LEFT (z<0, ITM
        % call); the negative shift (a<0, Ra=0) is accurate on the RIGHT (z>0,
        % OTM call). The blend weight w = sigmoid(5z) goes 0 -> 1 with z, so the
        % POS branch must carry weight (1-w) and the NEG branch weight w (same
        % convention as lewis_FFT_digital / ccdf_increment_FFT). Carrying w on
        % the pos branch (as before) inverted the blend and left kinks at z=+-1.
        price_grid              = zeros(size(z_grid));
        price_grid(z_grid < -1) = price_clean_pos(z_grid < -1);
        price_grid(z_grid >  1) = price_clean_neg(z_grid >  1);

        idx_blend             = (z_grid >= -1) & (z_grid <= 1);
        w                     = 1 ./ (1 + exp(-5 * z_grid(idx_blend)));
        price_grid(idx_blend) = (1 - w) .* price_clean_pos(idx_blend) + ...
                                      w  .* price_clean_neg(idx_blend);
    else
        price_grid = price_clean_neg;
    end

    % --- Interpolate at requested strikes ---
    price = interp1(z_grid, price_grid, strikes, 'pchip', 'extrap');

end