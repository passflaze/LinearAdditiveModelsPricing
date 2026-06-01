# Archive — funzioni ritirate (non sul path MATLAB)

Questa cartella **non** viene aggiunta al path da nessuno `run_*.m`, quindi i file
qui dentro non vengono mai eseguiti né possono fare shadowing. Sono conservati
per riferimento storico. Per riattivarne uno, riportalo nella sua cartella
originaria.

## Cluster 1 — engine MA vecchio (sostituito da `FA_simulation.m` + `Distributions/`)

| File | Originale | Sostituito da |
|------|-----------|---------------|
| `FA_simulationOLD.m` | `Simulation/Simulation_MA/` | `Simulation/Simulation_MA/FA_simulation.m` |
| `cf_FA_MA.m` | `Simulation/Simulation_MA/` | `Distributions/conditional_cf_MA_FA.m` (stessa formula, firma `(u, params, scale_factor)`) |
| `cf_IA_MA.m` | `Simulation/Simulation_MA/` | `Distributions/cf_MA_IA.m` |
| `lewis_fft_algorithm_digital.m` | `Simulation/Simulation_MA/` | `Simulation/lewis_FFT_digital.m` |
| `simulate_from_cdfPROVA.m` | `Simulation/Simulation_MA/` | `Simulation/simulate_from_cdf.m` |
| `exact_MA_cdf.m` | `Simulation/Simulation_MA/` | `Simulation/Simulation_MA/exact_ma_increment_cdf.m` |
| `tail_adjustment.m` | `Simulation/Simulation_MA/` | `Simulation/tail_adjustment.m` (copia IDENTICA, evita shadowing sul path) |

`cf_FA_MA`, `cf_IA_MA`, `lewis_fft_algorithm_digital` erano richiamati **solo** da
`FA_simulationOLD`; gli altri erano già orfani (0 riferimenti).

## Cluster 2 — approccio AB/GL vecchio (sostituito da `cf_increment_AB/GL` + `lewis_FFT_digital`)

| File | Originale | Note |
|------|-----------|------|
| `ccdf_increment_FFT.m` | `Simulation/` | nessun chiamante di codice; rimpiazzato da `lewis_FFT_digital` (vedi commento in `pricing_fwd_start_MC.m`) |
| `model_marginal_cf.m` | `Simulation/` | usato **solo** da `ccdf_increment_FFT` |
