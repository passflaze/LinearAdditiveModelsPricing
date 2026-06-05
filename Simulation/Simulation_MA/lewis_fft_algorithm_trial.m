function [cdf_clean, z_grid] = lewis_fft_algorithm_trial(M, dz, shift,clean)
%LEWIS_FFT_ALGORITHM_TRIAL  Sandbox CDF inversion via FFT (Gaussian CF test).
%   Inverts a characteristic function with the Lewis (2001) formula on a single
%   damping shift and optionally cleans the tails. Used as a standalone trial
%   harness; the CF is hard-wired to a standard Gaussian inside the routine.
%
%   INPUTS:
%       M     - FFT grid exponent, N = 2^M (e.g. 12 to 15)
%       dz    - step size of the frequency grid (z-domain)
%       shift - damping factor (alpha) shifting the integration contour
%       clean - (logical) if true, apply tail_adjustment post-processing
%
%   OUTPUTS:
%       cdf_clean - CDF vector (strictly monotonic and refined if clean = true)
%       z_grid    - corresponding spatial grid
%
%   See also: TAIL_ADJUSTMENT, LEWIS_FFT_DIGITAL

    % --- FFT grid setup ---
    N = 2^(M);
    dx = (2*pi) / (N * dz);         
    
    zn = (dz * (N-1)) / 2;
    z1 = -zn;
    z_grid = z1 : dz : zn;
    
    xn = (dx * (N-1)) / 2;
    x1 = -xn;
    x_grid = x1 : dx : xn;
    
    % --- FFT inversion of the CF ---
    prefactor = dx * exp(-1i * x1 * (z_grid));
    preprefactor = -exp(+shift*z_grid) / (2*pi);

    % Evaluate the CF at the shifted frequency grid.
    x_grid_shifted = x_grid + 1i*shift;
    cf_gauss = @(u) exp(-0.5.*(u.^2));
    fourier_function1 = cf_gauss(x_grid_shifted) ;
    fourier_function = fourier_function1./ (1i*(x_grid_shifted));
    
    j_minus_1 = 0:N-1;
    input_fft = fourier_function .* exp(-1i * z1 * dx * j_minus_1);
    
    fft_cdf = fft(input_fft);
    cdf_clean = real(preprefactor.* prefactor .* fft_cdf);


    % --- PLOT: RAW FFT CDF ---
    figure('Name', 'Raw Numerical CDF');
    plot(z_grid, cdf_clean, 'r-', 'LineWidth', 1.5);
    title('Raw CDF from FFT (Before Tail Adjustment)');
    xlabel('Spatial Grid (z)');
    ylabel('Probability');
    grid on;
    % -------------------------
    % --- Tail convergence diagnostic ---
    % Check that the raw CDF converges to 0 on the left and 1 on the right;
    % abs() guards against negative ringing at the boundaries.
    left_tail_error = abs(cdf_clean(1));
    right_tail_error = abs(1 - cdf_clean(end));
    
    max_tail_error = max(left_tail_error, right_tail_error);
    tolerance = 1e-4;

    fprintf('\n--- CDF Boundary Convergence Check ---\n');
    fprintf('Max boundary deviation: %.4e (Threshold: %.4e)\n', max_tail_error, tolerance);

    if max_tail_error < tolerance
        fprintf('Status: PASSED. The spatial grid bounds are wide enough.\n');
    else
        fprintf('[WARNING] Status: FAILED. The CDF does not converge near 0 and 1 at the grid boundaries.\n');
        fprintf('          Consider increasing M or modifying dz to widen the spatial grid limits.\n');
    end
    fprintf('--------------------------------------\n\n');

    % --- Post-processing: tail extrapolation and refinement ---
    % Clean numerical oscillations on the spatial grid (z_grid), refining the
    % resolution by a factor of 10.
    if clean
        [cdf_clean, z_grid] = tail_adjustment(z_grid, cdf_clean, 10);
    
        % --- PLOT: CLEANED CDF ---
        figure('Name', 'Cleaned Monotonic CDF');
        plot(z_grid, cdf_clean, 'b-', 'LineWidth', 2);
        title('Cleaned and Interpolated CDF');
        xlabel('Spatial Grid (x)');
        ylabel('Probability');
        ylim([-0.05, 1.05]); % Slight margins to visualize the 0 and 1 bounds
        grid on;
        % -------------------------
    end

end