function [X_st, x_grid, CDF_clean] = simulate_increments(zk, CDF_grid, N_sim, seed)

    % --- 1) TRONCAMENTO PER MONOTONIA -----------------------------------
    d_CDF = diff(CDF_grid);
    is_monotone = d_CDF > 0;

    [~, idx_max] = max(CDF_grid);

    idx0 = idx_max;
    while idx0 > 1 && is_monotone(idx0 - 1)
        idx0 = idx0 - 1;
    end

    idxK = idx_max;
    while idxK < length(CDF_grid) && is_monotone(idxK)
        idxK = idxK + 1;
    end

    x_grid    = zk(idx0:idxK);
    CDF_clean = CDF_grid(idx0:idxK);

    % --- 2) TRONCAMENTO PER RANGE [0,1] ---------------------------------
    % Dentro la regione monotona, trovo il sottoblocco dove CDF ∈ [0,1]
    in_range = CDF_clean >= 0 & CDF_clean <= 1;
    ids      = find(in_range);

    if isempty(ids)
        error('simulate_increments:emptyRange', ...
              'Nessun punto della CDF in [0,1] nella regione monotona.');
    end

    % Blocco contiguo (primo:ultimo indice valido)
    x_grid    = x_grid(ids(1):ids(end));
    CDF_clean = CDF_clean(ids(1):ids(end));

    % --- 3) CHECK CODE SUFFICIENTEMENTE SOTTILI -------------------------
    tol = 1e-4;
    if CDF_clean(1) > tol
        warning('simulate_increments:leftTail', ...
                'Coda sinistra troncata troppo presto: P(x_b)=%.3e > %.0e', ...
                CDF_clean(1), tol);
    end
    if 1 - CDF_clean(end) > tol
        warning('simulate_increments:rightTail', ...
                'Coda destra troncata troppo presto: 1-P(x_e)=%.3e > %.0e', ...
                1 - CDF_clean(end), tol);
    end

    fprintf('Troncata: P(x_b)=%.4e   1-P(x_e)=%.4e\n', ...
            CDF_clean(1), 1 - CDF_clean(end));

    % --- 4) FORZA UNICITÀ PER INVERTIBILITÀ -----------------------------
    [CDF_clean, uniq_idx] = unique(CDF_clean, 'stable');
    x_grid = x_grid(uniq_idx);

    % --- 5) SIMULAZIONE VIA INVERSIONE ----------------------------------
    rng(seed);
    U    = rand(N_sim, 1);
    U    = min(max(U, CDF_clean(1)), CDF_clean(end));
    X_st = interp1(CDF_clean, x_grid, U, 'spline');
end