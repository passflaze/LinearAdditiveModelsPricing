function [X_T1_grid, C_T1_grid] = call_GL_FFT(alpha_GL, beta_GL, t1, sigma_t1, t2, sigma_t2, B_T1_T2)
 
    M = 14; 
    dz = 0.005;
    N = 2^M;
    dx = 2 * pi / (N * dz); 
    
    z1 = -dz * (N - 1) / 2;
    x1 = -dx * (N - 1) / 2;
    j = 0:N-1;
    
    zk = z1 + dz * j;  
    xk = x1 + dx * j;  
    
    a = -1/2; 
    
    phi = @(u) cf_increment_GL(u, alpha_GL, beta_GL, sigma_t1, t1, sigma_t2, t2);
    
    int_eval = phi(xk + 1i*a) ./ (1i*xk - a).^2;
    fj = int_eval .* exp(-1i * z1 * dx .* j);
    f_hat_raw = exp(a * zk) .* dx .* exp(-1i * x1 * zk) .* fft(fj);
    
    % Prezzo non scontato
    C = real(f_hat_raw / (2 * pi));
    
    valid_idx = (zk >= -4) & (zk <= 4);
    zk_cut = zk(valid_idx);
    C_cut = C(valid_idx);
    
    % Ribaltamento coerente solo sulla zona sana
    X_T1_grid = flip(-zk_cut); 
    C_T1_grid = flip(C_cut) * B_T1_T2;
end
