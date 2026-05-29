function price_MC = price_forward_start_MC(alpha, beta, sigma_s, s, sigma_t, t, ...
                                     F0, B0t, strike, N_sim, seed)
%FORWARD_START_MC  Prezzo MC di una forward start option in modello GL.
%
%   price_MC = FORWARD_START_MC(alpha, beta, sigma_s, s, sigma_t, t, ...
%                               F0, B0t, strike, N_sim)
%
%   Payoff: [S_t - strike * F(s,t)]^+
%   Decomposizione: f_{s,t} = X1,  f_{t,t} = X1 + X2
%
%   INPUT:
%     alpha, beta - parametri GL
%     sigma_s, s  - sigma e tempo della reset date
%     sigma_t, t  - sigma e tempo della expiry
%     F0          - F(0, t), forward iniziale con expiry t
%     B0t         - discount factor B(0, t)
%     strike      - K2
%     N_sim       - numero di simulazioni
%
%   OUTPUT:
%     price_MC    - prezzo Monte Carlo

    % Incremento 0 -> s
    [CDF1, zk1] = lewis_fft_cdf(@cf_increment_GL, alpha, beta, ...
                                0, 0, sigma_s, s);
    [X1, ~, ~]  = simulate_increments(zk1, CDF1, N_sim, seed);

    % Incremento s -> t
    [CDF2, zk2] = lewis_fft_cdf(@cf_increment_GL, alpha, beta, ...
                                sigma_s, s, sigma_t, t);
    [X2, ~, ~]  = simulate_increments(zk2, CDF2, N_sim, seed);

    % Ricostruzione forward
    F_st = F0 + X1;          % F(s, t) %forse stiamo calcolando fss con in questo modo
    S_t  = F0 + X1 + X2;     % S_t = F(t, t)

    % Prezzo MC
    payoff   = max(S_t - strike .* F_st, 0);
    price_MC = B0t * mean(payoff);
end