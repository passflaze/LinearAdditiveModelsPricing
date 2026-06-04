function [cdf_grid, z_grid] = lewis_FFT_digital(cf, M, dz, params, scale_factors, clean, model, doubleshift, doplot)
% LEWIS_FFT_ALGORITHM_DIGITAL  Computes digital put prices via FFT inversion
%                               of the CF using the Lewis (2001) formula.
%
%   A sigmoid blend joins the two shifted CDF grids in the transition
%   region [-1, 1]. An optional tail adjustment enforces monotonicity.
%
% INPUTS:
%   cf          - (function handle) CF of the increment: cf(u, params, scale_factor)
%   M           - (scalar) grid size exponent, N = 2^M (typical: 12-15)
%   dz          - (scalar) step size in the z-domain (log-strike grid)
%   params      - (vector) model parameters passed to cf
%   scale_factor- [scale_t1, scale_t2] scaling factors at t1 and t2
%   z_grid_std  - (vector) standardised z-grid for spline resampling
%   clean       - (logical) if true, applies tail_adjustment post-processing
%   model       - (string) model identifier: 'MA', 'GL', or 'AB'
%   doubleshift - (logical) if true, uses both shifts and blends the grids
%
% OUTPUTS:
%   cdf_grid    - (vector) digital put price (CDF) grid, optionally cleaned
%   z_grid      - (vector) corresponding log-strike spatial grid


    % Plotting + verbose boundary report are OFF by default: this routine is
    % called repeatedly (MC, every model/maturity), so unconditional figures
    % and prints would spam. Pass doplot = true for the diagnostic figures.
    if nargin < 9 || isempty(doplot)
        doplot = false;
    end

    % Use the largest scale ONLY for the damping shift. 
    scale_factor = max(scale_factors);
    % --- FFT grid setup ---
    N  = 2^M;
    dx = (2*pi) / (N * dz);

    zn     = (dz * (N-1)) / 2;
    z1     = -zn;
    z_grid = z1 : dz : zn;

    xn     = (dx * (N-1)) / 2;
    x1     = -xn;
    x_grid = x1 : dx : xn;

    % --- Damping shifts ---
    switch model
        case 'MA'
            % params(1) = alpha, params(2) = beta
            shift_pos =  0.5 * params(1) / scale_factor;
            shift_neg = -0.5 * params(2) / scale_factor;

        case 'GL'
            % params(1) = alpha, params(2) = beta
            shift_pos =  0.5 * params(1) / scale_factor;
            shift_neg = -0.5 * params(2) / scale_factor;

        case 'AB'
            % params(1) = k (kappa), params(2) = eta
            % p_t^+/- = (1/scale_t2) * ( +-eta + sqrt(eta^2 + 1/kappa) )
            % (alpha = 0.5 => 2*(1-alpha)/kappa = 1/kappa)
            kappa = params(1);
            eta   = params(2);
            shift_pos =  0.5 * (-eta + sqrt(eta^2 + 1/kappa)) / scale_factor;
            shift_neg = -0.5 * (eta + sqrt(eta^2 + 1/kappa)) / scale_factor;
    end

    % --- Positive shift: FFT inversion ---
    j_minus_1        = 0:N-1;
    prefactor        = dx * exp(-1i * x1 * z_grid);
    preprefactor_pos = -exp(shift_pos * z_grid) / (2*pi);

    x_grid_shifted_pos   = x_grid + 1i * shift_pos;
    fourier_function_pos = cf(x_grid_shifted_pos, params, scale_factors) ./ ...
                           (1i * x_grid_shifted_pos);
    input_fft_pos        = fourier_function_pos .* exp(-1i * z1 * dx * j_minus_1);

    fft_cdf_pos   = real(prefactor .* fft(input_fft_pos));
    % fft_cdf_pos   = interp1(z_grid, fft_cdf_pos, z_grid_std, 'spline');
    % preprefactor_pos = interp1(z_grid, preprefactor_pos, z_grid_std, 'spline');
    cdf_clean_pos = preprefactor_pos .* fft_cdf_pos;

    % --- Negative shift + blending (doubleshift mode) ---
    if doubleshift
        preprefactor_neg = -exp(shift_neg * z_grid) / (2*pi);

        x_grid_shifted_neg   = x_grid + 1i * shift_neg;
        fourier_function_neg = cf(x_grid_shifted_neg, params, scale_factors) ./ ...
                               (1i * x_grid_shifted_neg);
        input_fft_neg        = fourier_function_neg .* exp(-1i * z1 * dx * j_minus_1);

        fft_cdf_neg      = real(prefactor .* fft(input_fft_neg));
        % fft_cdf_neg      = interp1(z_grid, fft_cdf_neg,      z_grid_std, 'spline');
        % preprefactor_neg = interp1(z_grid, preprefactor_neg, z_grid_std, 'spline');
        cdf_clean_neg    = 1 + preprefactor_neg .* fft_cdf_neg;
        z_grid_std = z_grid;

        % Each one-shift reconstruction is a CDF (digital-put price) and must
        % lie in [0,1]. At the extreme grid edges (|z| ~ N*dz/2) the damping
        % prefactor exp(shift*z) of the OPPOSITE side overflows, so the raw
        % grid blows up there (e.g. ~1e23 for the gamma-based GL CF). Clamp
        % BEFORE the blend: this is a no-op in the valid central region but
        % removes the overflow, so the tanh blend (and the boundary-convergence
        % diagnostic below) stay finite. tail_adjustment then splices the tails.
        cdf_clean_pos = min(max(cdf_clean_pos, 0), 1);
        cdf_clean_neg = min(max(cdf_clean_neg, 0), 1);

        % Tanh blend
        scale    = (z_grid_std(end) - z_grid_std(1)) / 20;
        w        = 0.5 * (1 + tanh(z_grid_std / scale));
        cdf_grid = (1 - w) .* cdf_clean_pos + w .* cdf_clean_neg;
    else
        cdf_grid = min(max(cdf_clean_pos, 0), 1);
        z_grid_std = z_grid;
    end

    % --- PLOT: RAW FFT CDF WITH DYNAMIC TITLE (optional) ---
    if doplot
        figure('Name', ['Raw Numerical CDF ' model]);
        plot(z_grid_std, cdf_grid, 'r-', 'LineWidth', 1.5);
        title(['Raw CDF from FFT (Before Tail Adjustment) ' model]);
        xlabel('Spatial Grid (z)');
        ylabel('Probability');
        grid on;
    end

    % --- Tail convergence diagnostic: warn only if the bounds are too tight ---
    left_tail_error  = abs(cdf_grid(1));
    right_tail_error = abs(1 - cdf_grid(end));
    max_tail_error   = max(left_tail_error, right_tail_error);
    tolerance        = 1e-4;
    if max_tail_error >= tolerance
        warning('lewis_FFT_digital:BoundaryConvergence', ...
            ['%s: CDF boundary deviation %.2e exceeds %.0e. ', ...
             'Consider increasing M or modifying dz.'], ...
            model, max_tail_error, tolerance);
    end

    % --- Optional tail adjustment ---
     if clean
        [cdf_grid, z_grid] = tail_adjustment(z_grid_std, cdf_grid, 10);

        if doplot
            figure('Name', ['Cleaned Monotonic CDF ' model]);
            plot(z_grid, cdf_grid, 'b-', 'LineWidth', 2);
            title(['Cleaned and Interpolated CDF ' model]);
            xlabel('Spatial Grid (x)');
            ylabel('Probability');
            ylim([-0.05, 1.05]);
            grid on;
        end
    end


end