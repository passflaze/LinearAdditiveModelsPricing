# Linear Additive Models for Option Pricing

A MATLAB quantitative library for calibrating, simulating, pricing, and hedging
options on **WTI crude-oil futures** using three **Linear Additive Models**:

| Model | Increment law | Key reference |
|-------|---------------|---------------|
| **AB** — Additive Bachelier | Normal-Inverse-Gaussian (tempered-stable, α = ½) | Baviera & Massaria (2026) |
| **MA** — Minimal Additive | Asymmetric Laplace (double-exponential) | Carr & Torricelli; Baviera et al. |
| **GL** — Generalized Logistic | Generalized-logistic / additive logistic | Carr & Torricelli (2021) |

The library covers the full lifecycle of a structured-products desk: it bootstraps
the discount and forward curves from listed options, fits each additive model to the
implied-volatility surface, prices vanilla and **exotic path-dependent** options
(Call-on-Call, Put-on-Put, Chooser, forward-start) via both **closed-form** and
**Lewis-FFT Monte-Carlo** engines, and finally builds and back-tests **dynamic
Delta-Gamma-Vega hedges** with realistic transaction costs.

> **Reference snapshot:** all results are computed on the WTI option surface as of
> the value date **2 June 2020**, across 9 listed maturities (JUN20 → DEC22).

---

## Why "Linear Additive" models?

Classical Lévy models assume *stationary* increments: the law of the price change
over `[t, t+Δ]` depends only on `Δ`. Real option surfaces violate this — the
implied-volatility smile flattens with maturity in a way stationary models cannot
reproduce with a single parameter set.

**Additive** processes relax stationarity: increments stay independent but their law
is allowed to change over time. In this library that time dependence is carried by a
single **scale factor**

```
scale_factor(t) = sigma_t * sqrt(t)
```

so each model is described by just **two shape parameters** (plus the ATM volatility
term structure), yet still fits the whole surface. The three models differ only in
the *shape* of the standardized increment — Gaussian-like core with semi-heavy tails
(AB), exponential tails (MA), or logistic tails (GL) — which is exactly the
comparison drawn in `run_ex0`.

---

## Repository layout

```
LinearAdditiveModelsPricing/
├── run_project2A.m          ← MAIN entry point: orchestrates Ex0 → Ex6 below
├── ex2.m                    Calibration-only convenience wrapper
│
├── Run_Functions/           High-level stage drivers, called by the main script
│   ├── run_ex0.m            Ex0  – analytical PDF comparison of the 3 models
│   ├── run_ex2.m            Ex1+2 – curve bootstrap + surface calibration
│   ├── run_ex3.m            Ex3  – forward-start pricing (Analytic vs MC vs FFT)
│   ├── run_ex4.m            Ex4  – exotic pricing (CoC / PoP / Chooser)
│   └── run_ex6.m            Ex6  – dynamic Delta-Gamma-Vega hedging backtest
│
├── Utilities/               Curve bootstrap, ATM vol, Bachelier IV, data readers
├── Distributions/           Characteristic functions, PDFs and CDFs per model
├── Calibration/             Surface-fitting objective functions and I0 normalizers
│   ├── Calibration_AB/      AB-specific pricing & calibration helpers
│   ├── Calibration_MA/      MA-specific pricing & calibration helpers
│   └── Calibration_GL/      GL-specific pricing & calibration helpers
├── Simulation/              Lewis-FFT engine + CDF-inversion Monte-Carlo sampler
│   └── Simulation_MA/       Exact MA increment sampler & moment diagnostics
├── Pricing/
│   ├── Analytic/            Closed-form vanilla, forward-start & exotic prices
│   └── MC/                  Monte-Carlo pricers for the same payoffs
├── Hedging/                 Greeks, hedge construction, costs, and backtester
│
├── Data/                    WTI call/put price grids + futures expiry calendar
├── Biblio/                  Reference papers (the project's source of truth)
└── Archive/                 Retired functions kept for reference (NOT on the path)
```

> **Convention:** every vector — strikes, maturities, prices, parameters — is a
> **column vector** throughout the codebase, and all element-wise vs. matrix
> operators are chosen accordingly.

---

## The pipeline

Running [run_project2A.m](run_project2A.m) executes six stages in sequence.
Calibrated parameters flow forward through the shared `params` and `market` structs —
**nothing is hard-coded downstream**.

### Ex 0 — Model PDF comparison · [run_ex0.m](Run_Functions/run_ex0.m)
Plots the analytical densities of the three standardized increments (Asymmetric
Laplace, NIG, Generalized Logistic) on linear, semi-log, and extreme-tail scales,
benchmarked against the Normal. An optional 1M-path Monte-Carlo overlay validates the
AB (NIG) increment generator against its closed form.

### Ex 1 & 2 — Curve bootstrap + calibration · [run_ex2.m](Run_Functions/run_ex2.m)
1. **Bootstrap** ([bootstrap.m](Utilities/bootstrap.m)): for each maturity, the
   discount factor `D(t,T)` and forward `F(t,T)` are extracted by OLS regression on
   put-call parity `C - P = D·(F - K)` (synthetic-forward technique of Azzone &
   Baviera). `R²` reports how tightly parity holds.
2. **ATM vol term structure**: the at-the-money price is interpolated and inverted to
   a **Bachelier** (normal) implied volatility.
3. **OTM surface filter**: only out-of-the-money options in the dollar-moneyness band
   `x = K − F ∈ [−30$, 30$]` are kept; OTM puts are converted to synthetic calls via
   parity (per Baviera & Massaria, Sec. 4.1).
4. **Calibration**: each model's two parameters are fit to the filtered call surface
   by least squares — AB and GL via `fmincon`, MA via 1-D `fminbnd` with a fixed
   scale gauge. Returns:
   - `params.AB = [k; eta]`, `params.MA = [alpha; beta]`, `params.GL = [alpha; beta]`
   - `market` — all curves, vols, year-fractions and the filtered surface.

### Ex 3 — Forward-start pricing · [run_ex3.m](Run_Functions/run_ex3.m)
Prices a forward-start option struck at `K2` over the window `[t1, t2]`
(here T1 = index 2, T2 = index 4). For each model it cross-checks **three independent
engines**: the closed-form price, a CDF-inversion **Monte Carlo**, and a uniform
**Lewis-FFT** inversion — agreement is reported in basis points. For MA it also
verifies the analytic price against direct numerical integration of the exact
increment CDF.

### Ex 4 — Exotic pricing · [run_ex4.m](Run_Functions/run_ex4.m)
Prices three compound / path-dependent payoffs — **Call-on-Call (CoC)**,
**Put-on-Put (PoP)**, and **Chooser** — by Monte Carlo, with closed-form prices where
available (fully analytic for MA). It ships two built-in **sanity checks**:
- `K1 = 0` collapses a Call-on-Call to the underlying vanilla ATM call;
- when `K2 = F(t0,T2)`, the Chooser equals the sum of two ATM calls (symmetry).

An optional **smart-extrapolation** test quantifies the pricing error introduced by
truncating the simulation CDF, comparing analytic-tail reconstruction against a brute
truncation.

### Ex 6 — Dynamic hedging · [run_ex6.m](Run_Functions/run_ex6.m)
Builds a static hedge that neutralizes a chosen set of Greeks (any subset of
**Delta / Gamma / Vega**) of an exotic portfolio, using a user-defined basket of
vanilla calls, puts and futures. Greeks are computed by finite differences (bump the
forward and the vol term structure); hedge weights solve the linear Greek-matching
system ([build_hedge.m](Hedging/build_hedge.m)); and opening **transaction costs**
are charged as a bid-ask spread (1 bp futures, 4 bp options). The position is then
**back-tested** over successive rebalancing dates (Tuesdays), tracking P&L, residual
Greeks and cumulative cost. The main script runs six scenarios (long/short Chooser,
CoC+PoP baskets, Delta- / Gamma- / Vega-targeted hedges).

---

## Numerical engines

- **Lewis-FFT** ([lewis_FFT_call.m](Simulation/lewis_FFT_call.m),
  [lewis_FFT_digital.m](Simulation/lewis_FFT_digital.m)) — call and digital prices by
  Fourier inversion of the characteristic function, using **dual damping shifts**
  blended across a transition region to keep the inversion stable in both wings.
- **CDF-inversion Monte Carlo** ([simulate_from_cdf.m](Simulation/simulate_from_cdf.m),
  [smart_cdf_extrapolation.m](Simulation/smart_cdf_extrapolation.m)) — samples
  increments by inverting an FFT-built CDF, with analytic-tail extrapolation to avoid
  losing mass at truncation.
- **Forward-start CF rescaling** — the increment CF of the `T2`-forward is built as
  `phi_t2(u) / phi_t1(fwd_factor · u)` with `fwd_factor = B(0,T1)/B(0,T2)`
  (Lemma 2, `Biblio/Forward.pdf`).
- **Exact MA sampler** ([FA_simulation.m](Simulation/Simulation_MA/FA_simulation.m))
  — analytic isolation of the discrete point-mass to suppress Gibbs oscillations.

---

## Getting started

**Requirements:** MATLAB R2020a+ with the Optimization, Statistics & Machine
Learning, and Financial toolboxes (`fmincon`, `fminbnd`, `ksdensity`, `yearfrac`).

**Run the full pipeline:**

```matlab
>> run_project2A
```

The script adds all sub-folders to the path, runs Ex0 → Ex6, and prints calibration
reports, pricing tables (in bps), and hedging diagnostics to the console. To see the
figures, flip the per-stage flags at the top of each block:

```matlab
plot_ex0     = true;     % Ex0 density plots
opts_ex2.plot    = true; % bootstrap curves, smile, term-structure
opts_ex2.verbose = true; % full diagnostic report
opts_ex3.plot    = true; % forward-start price curves
opts_ex4.plot    = true; % CoC / PoP price vs strike
```

**Run a single stage** — each `run_exN` is self-contained and will calibrate
internally if no `params`/`market` are passed:

```matlab
[params, market]   = run_ex2(struct('verbose', true, 'plot', true));
LA_results_es4     = run_ex4(params, market, opts_ex4);
```

The calibration-only wrapper [ex2.m](ex2.m) is a convenient shortcut to inspect the
fitted parameters and implied distributions.

---

## Bibliography (`Biblio/`)

This project implements **specific, non-standard formulations**; the papers in
[Biblio/](Biblio/) are the authoritative source for every formula, characteristic
function, and parameter constraint:

1. Azzone & Baviera — *A fast Monte Carlo scheme for additive processes*
2. Azzone & Baviera — *Synthetic forwards and cost of funding in the equity
   derivative market* (curve bootstrap)
3. Baviera & Massaria (2026) — *Additive Bachelier* (surface filtering, AB model)
4. Carr & Torricelli — *Additive logistic processes in option pricing* (GL model)
5. Baviera & Manzoni (2026) — *FGMC* (fast generalized Monte Carlo)
6. `Forward.pdf` — forward-start CF rescaling (Lemma 2)
