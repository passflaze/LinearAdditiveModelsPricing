function increments = FA_simulation(N_sim, M, dz, delta_mu, pt_plus, pt_minus, ps_plus, ps_minus,spline,activity,doubleshift, params, scale_factor, doPlot)
%FA_SIMULATION Simulate a Minimal Additive (MA) increment by CDF inversion.
%   For the finite-activity case it uses a Poisson mixture: a Bernoulli trial
%   selects which paths jump, and the jump magnitudes are drawn (only for the
%   active paths) from the FFT-reconstructed CDF. For the infinite-activity
%   case every path is sampled directly from the inverted CDF.
%
%   INPUTS:
%       N_sim        - number of Monte Carlo paths to generate
%       M            - FFT grid exponent (N = 2^M)
%       dz           - step size of the spatial grid
%       delta_mu     - deterministic drift (martingale correction) for [s, t]
%       pt_plus      - right-tail decay parameter at terminal time t
%       pt_minus     - left-tail decay parameter at terminal time t
%       ps_plus      - right-tail decay parameter at initial time s
%       ps_minus     - left-tail decay parameter at initial time s
%       spline       - (logical) PCHIP (true) vs linear inversion in simulate_from_cdf
%       activity     - (string) 'finite' or 'infinite'
%       doubleshift  - (logical) dual-shift FFT inversion in lewis_FFT_digital
%       params       - [alpha; beta] MA parameters passed to the CF
%       scale_factor - [scale_t1; scale_t2] scaling factors at t1 and t2
%       doPlot       - (optional, default false) plot the reconstructed CDF
%
%   OUTPUT:
%       increments - [N_sim x 1] vector of simulated process increments
%
%   See also: POINTMASSCALCULATOR, SIMULATE_FROM_CDF, LEWIS_FFT_DIGITAL

   if nargin < 14 || isempty(doPlot)
    doPlot = false;
    end

   if strcmp(activity, 'finite')
        % Jump probability from the total intensity: P(K>0) = 1 - exp(-lambda).
        lambda = pointmasscalculator(pt_plus, pt_minus, ps_plus, ps_minus);
        p_jump = 1 - exp(-lambda);

        % Bernoulli split: which paths jump.
        jump_mask = rand(N_sim, 1) < p_jump;
        N_jumps = sum(jump_mask);

        jump_component = zeros(N_sim, 1);

        % Invert the CDF and draw magnitudes only for the paths that jumped.
        if N_jumps > 0
            [cdf_clean, x_grid_fine] = lewis_FFT_digital(@conditional_cf_MA_FA,M, dz, params,scale_factor, 1,'MA',doubleshift, doPlot);
            simulated_jumps = simulate_from_cdf(cdf_clean, x_grid_fine,spline,N_jumps);
            jump_component(jump_mask) = simulated_jumps;
        end
        increments = jump_component + delta_mu;
   else
        [cdf_clean, x_grid_fine] = lewis_FFT_digital(@cf_MA_IA,M, dz, params,scale_factor, 1,'MA',doubleshift, doPlot);
        increments = simulate_from_cdf(cdf_clean, x_grid_fine,spline,N_sim);
        
   end
 
end

