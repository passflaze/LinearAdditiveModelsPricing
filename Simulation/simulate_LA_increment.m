function [Z, x_fine, cdf_fine] = simulate_LA_increment(spec, T1, T2, ...
                                          sigma_T1, sigma_T2, fwd_factor, N_sim)
%SIMULATE_LA_INCREMENT  Draw N_sim samples of a Linear Additive increment
%   (AB or GL) via the shared pipeline:
%       ccdf_LA_FFT  ->  tail_adjustment  ->  simulate_from_cdf
%   i.e. the two routines the project already validated (tail_adjustment and
%   simulate_from_cdf) are reused unchanged; only the raw CDF generator is the
%   unified ccdf_LA_FFT.
%
%   [Z, X_FINE, CDF_FINE] = SIMULATE_LA_INCREMENT(SPEC, T1, T2, SIGMA_T1, ...
%                                                 SIGMA_T2, FWD_FACTOR, N_SIM)
%
%   For the MARGINAL increment 0 -> T pass T1 = 0, sigma_T1 = 0, fwd_factor = 1
%   (then T2/sigma_T2 are the target maturity/scale).
%   For the CONDITIONAL increment T1 -> T2 pass the Lemma-2 fwd_factor.
%
%   INPUTS
%     spec       : model spec struct (la_model_spec)
%     T1, T2     : reset and maturity year fractions (T1 = 0 -> marginal)
%     sigma_T1   : model scale at T1 (0 if marginal)
%     sigma_T2   : model scale at T2
%     fwd_factor : Lemma-2 forward rescaling B(0,T1)/B(0,T2) (1 if marginal)
%     N_sim      : number of samples
%
%   OUTPUTS
%     Z        : N_sim x 1 simulated increments
%     x_fine   : refined spatial grid of the cleaned CDF (column)
%     cdf_fine : cleaned, strictly monotone CDF on x_fine (column) -- reusable
%                for an FFT/integration reference price.

    if nargin < 6 || isempty(fwd_factor)
        fwd_factor = 1;
    end

    % --- auto-size the spatial grid to the increment std ----------------------
    % Bachelier scaling: dispersion = sigma * sqrt(T). The standardized variance
    % spec.varStd converts it to the increment variance. 10 std pushes each tail
    % mass below ~1e-4 (the rule-of-thumb threshold used in tail diagnostics).
    sT2 = sigma_T2 * sqrt(T2);
    if T1 == 0
        var_inc = spec.varStd * sT2^2;
    else
        sT1     = sigma_T1 * sqrt(T1);
        var_inc = spec.varStd * (sT2^2 - fwd_factor^2 * sT1^2);
    end
    std_inc = sqrt(max(var_inc, eps));

    n_grid = 2000;
    x_grid = linspace(-10*std_inc, 10*std_inc, n_grid)';

    % --- raw conditional CDF via unified Lewis-FFT ---------------------------
    cdf_raw = ccdf_LA_FFT(spec, T1, T2, sigma_T1, sigma_T2, x_grid, fwd_factor);

    % --- clean / refine / exponential-tail splice (validated routine) --------
    [cdf_fine, x_fine] = tail_adjustment(x_grid, cdf_raw, 10);
    cdf_fine = cdf_fine(:);
    x_fine   = x_fine(:);

    % --- inverse-CDF sampling (validated routine, spline/pchip branch) -------
    Z = simulate_from_cdf(cdf_fine, x_fine, true, N_sim);
    Z = Z(:);
end
