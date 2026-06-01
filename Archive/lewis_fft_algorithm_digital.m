function [cdf_grid, z_grid] = lewis_fft_algorithm_digital(cf,M, dz, shift_pos, shift_neg, pt_plus, pt_minus, ps_plus, ps_minus, clean, activity, doubleshift)
%LEWIS_FFT_ALGORITHM Computes the CDF of a Minimal Additive process/increment via FFT.
%   [CDF_GRID, Z_GRID] = LEWIS_FFT_ALGORITHM(M, DZ, SHIFT_POS, SHIFT_NEG, 
%   PT_PLUS, PT_MINUS, PS_PLUS, PS_MINUS, CLEAN, ACTIVITY) uses the Fast Fourier 
%   Transform to invert the characteristic function of an MA process. It applies 
%   dual damping shifts to ensure integrability and can use a post-processing 
%   tail adjustment to guarantee a strict, non-oscillatory CDF profile.
%
%   Inputs:
%       M         - Controls grid size, where N = 2^M (e.g., 12 to 15)
%       dz        - Step size for the frequency grid (z-domain)
%       shift_pos - Positive damping factor for the right-hand integration
%       shift_neg - Negative damping factor for the left-hand integration
%       pt_plus   - Right-tail decay parameter at terminal time t
%       pt_minus  - Left-tail decay parameter at terminal time t
%       ps_plus   - Right-tail decay parameter at initial time s
%       ps_minus  - Left-tail decay parameter at initial time s
%       clean     - Boolean flag (true/false) to activate tail_adjustment
%       activity  - String identifier: 'finite' or 'infinite'
%       doubleshift- Boolean flag whether to apply or not the double shift 
%
%   Outputs:
%       cdf_grid  - The evaluated (and optionally cleaned) CDF vector
%       z_grid    - The corresponding high-resolution spatial grid

    % =========================================================================
    % STEP 1: DYNAMIC TITLE STRING FOR GRAPHICS
    % =========================================================================
    if strcmp(activity, 'finite')
        activity_title = ' (Finite Activity)';
    else
        activity_title = ' (Infinite Activity)';
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
    preprefactor_pos = -exp(+shift_pos*z_grid) / (2*pi);
    
    % Evaluate the conditional CF at the shifted frequency grid
    x_grid_shifted_pos = x_grid + 1i*shift_pos;
    
 
    fourier_function_pos = cf(x_grid_shifted_pos, pt_plus, pt_minus, ps_plus, ps_minus) ./ (1i*(x_grid_shifted_pos));


    j_minus_1 = 0:N-1;
    input_fft_pos = fourier_function_pos .* exp(-1i * z1 * dx * j_minus_1);
    
    fft_cdf_pos = real(prefactor .* fft(input_fft_pos));
    cdf_clean_pos = preprefactor_pos .* fft_cdf_pos;

    if doubleshift 
        preprefactor_neg = -exp(+shift_neg*z_grid) / (2*pi);
        x_grid_shifted_neg = x_grid + 1i*shift_neg;
        fourier_function_neg = cf(x_grid_shifted_neg, pt_plus, pt_minus, ps_plus, ps_minus) ./ (1i*(x_grid_shifted_neg));
        input_fft_neg = fourier_function_neg .* exp(-1i * z1 * dx * j_minus_1);
        fft_cdf_neg = real(prefactor .* fft(input_fft_neg));
        cdf_clean_neg = 1 + preprefactor_neg .* fft_cdf_neg;

        % Union logic
        cdf_grid = zeros(size(z_grid));
        idx_sx = z_grid < -1;
        cdf_grid(idx_sx) = cdf_clean_pos(idx_sx);
        idx_dx = z_grid > 1;
        cdf_grid(idx_dx) = cdf_clean_neg(idx_dx);
        idx_blend = (z_grid >= -1) & (z_grid <= 1);
        c_sx = cdf_clean_pos(idx_blend);
        c_dx = cdf_clean_neg(idx_blend);
        exp_term = 1 ./ (1 + exp(-5 * z_grid(idx_blend)));
        cdf_grid(idx_blend) = exp_term .* c_sx + (1 - exp_term) .* c_dx;
    else
        cdf_grid = cdf_clean_pos;
    end
        
    % --- PLOT: RAW FFT CDF WITH DYNAMIC TITLE ---
    figure('Name', ['Raw Numerical CDF' activity_title]);
    plot(z_grid, cdf_grid, 'r-', 'LineWidth', 1.5);
    title(['Raw CDF from FFT (Before Tail Adjustment)' activity_title]);
    xlabel('Spatial Grid (z)');
    ylabel('Probability');
    grid on;
    
    % =========================================================================
    % STEP 3.5: TAIL CONVERGENCE DIAGNOSTIC CHECK
    % =========================================================================
    left_tail_error = abs(cdf_grid(1));
    right_tail_error = abs(1 - cdf_grid(end));
    
    max_tail_error = max(left_tail_error, right_tail_error);
    tolerance = 1e-4;
    fprintf('\n--- CDF Boundary Convergence Check%s ---\n', activity_title);
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
    if clean
        [cdf_grid, z_grid] = tail_adjustment(z_grid, cdf_grid, 10);
    
        % --- PLOT: CLEANED CDF WITH DYNAMIC TITLE ---
        figure('Name', ['Cleaned Monotonic CDF' activity_title]);
        plot(z_grid, cdf_grid, 'b-', 'LineWidth', 2);
        title(['Cleaned and Interpolated CDF' activity_title]);
        xlabel('Spatial Grid (x)');
        ylabel('Probability');
        ylim([-0.05, 1.05]); 
        grid on;
    end
end