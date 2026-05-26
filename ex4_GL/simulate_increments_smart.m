function [X_sim, x_grid_clean, CDF_clean] = simulate_increments_smart(x_grid, CDF_grid, N_sim, seed, p_minus, p_plus)

    rng(seed);
    U = rand(N_sim, 1);
    
    % 1. Pulizia della griglia FFT (togliere i NaN e i valori non strettamente crescenti)
    valid_idx = ~isnan(CDF_grid) & ~isinf(CDF_grid) & (CDF_grid > 0) & (CDF_grid < 1);
    x_valid = x_grid(valid_idx);
    CDF_valid = CDF_grid(valid_idx);
    
    % Assicuriamoci che la CDF sia strettamente crescente (necessario per interp1)
    [CDF_clean, unique_idx] = unique(CDF_valid, 'stable');
    x_grid_clean = x_valid(unique_idx);
    
    % 2. Identificazione dei bordi
    x_min = x_grid_clean(1);
    F_min = CDF_clean(1);
    
    x_max = x_grid_clean(end);
    F_max = CDF_clean(end);
    
    % 3. Allocazione del vettore dei risultati
    X_sim = zeros(N_sim, 1);
    
    % --- SMART EXTRAPOLATION ---
    
    % A. Coda Sinistra
    idx_left = U < F_min;
    X_sim(idx_left) = x_min + (1 / p_minus) * log(U(idx_left) / F_min);
    
    % B. Coda Destra
    idx_right = U > F_max;
    X_sim(idx_right) = x_max - (1 / p_plus) * log((1 - U(idx_right)) / (1 - F_max));
    
    % C. Zona Centrale (FFT)
    idx_center = (U >= F_min) & (U <= F_max);
    X_sim(idx_center) = interp1(CDF_clean, x_grid_clean, U(idx_center), 'pchip'); % pchip evita oscillazioni
end