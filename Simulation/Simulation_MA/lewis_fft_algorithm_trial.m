function [cdf_clean, z_grid] = lewis_fft_algorithm_trial(M, dz, shift,clean)
%INVERT_MA_CDF_FFT Computes the CDF of a Minimal Additive increment via FFT.
%   [CDF_CLEAN, Z_GRID_FINE] = INVERT_MA_CDF_FFT(M, DZ, SHIFT, PT_PLUS, PT_MINUS, 
%   PS_PLUS, PS_MINUS) uses the Fast Fourier Transform to invert the 
%   characteristic function of a finite-activity MA process. It applies a 
%   damping shift to ensure integrability and uses a post-processing 
%   spline/exponential tail adjustment to guarantee a strict CDF profile.
%
%   Inputs:
%       M         - Controls grid size, where N = 2^M (e.g., 12 to 15)
%       dz        - Step size for the frequency grid (z-domain)
%       shift     - Damping factor (alpha) to shift the integration contour
%       pt_plus   - Right-tail decay parameter at terminal time t
%       pt_minus  - Left-tail decay parameter at terminal time t
%       ps_plus   - Right-tail decay parameter at initial time s
%       ps_minus  - Left-tail decay parameter at initial time s
%
%   Outputs:
%       cdf_clean   - The strictly monotonic, cleaned CDF vector
%       z_grid_fine - The corresponding high-resolution spatial grid
%
%   See also: CF_FA_MA, TAIL_ADJUSTMENT

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
    preprefactor = -exp(+shift*z_grid) / (2*pi);
    
    % Evaluate the conditional CF at the shifted frequency grid
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
% =========================================================================
    % STEP 3.5: TAIL CONVERGENCE DIAGNOSTIC CHECK
    % =========================================================================
    % Check if the raw CDF properly converges to 0 on the left and 1 on the right.
    % We use abs() to account for potential negative ringing at the boundaries.
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

    % =========================================================================
    % STEP 4: POST-PROCESSING (TAIL EXTRAPOLATION & SPLINE INTERPOLATION)
    % =========================================================================
    % Clean numerical oscillations using the spatial grid (z_grid), NOT the 
    % frequency grid (x_grid), refining the resolution by a factor of 10.
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