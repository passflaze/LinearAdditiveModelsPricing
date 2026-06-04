function increments = FA_simulation(N_sim, M, dz, delta_mu, pt_plus, pt_minus, ps_plus, ps_minus,spline,activity,doubleshift, params, scale_factor, doPlot)
%SIMULATE_FA_INCREMENT Command center for simulating Finite Activity MA jumps.
%   INCREMENTS = SIMULATE_FA_INCREMENT(N_SIM, M, DX, SHIFT, DELTA_MU, ...)
%   orchestrates the full simulation of an increment for a Minimal Additive 
%   process. It uses a Poisson mixture approach: it first determines which 
%   paths experience at least one jump via a Bernoulli trial, and then 
%   generates the jump magnitudes exclusively for those active paths using 
%   the FFT-based CDF inversion.
%
%   Inputs:
%       N_sim     - Number of Monte Carlo paths to generate
%       M         - Controls grid size for FFT (N = 2^M)
%       dx        - Step size for the spatial grid
%       shift     - Damping factor for FFT contour integration
%       delta_mu  - The deterministic drift (martingale correction) for [s, t]
%       pt_plus   - Right-tail decay parameter at terminal time t
%       pt_minus  - Left-tail decay parameter at terminal time t
%       ps_plus   - Right-tail decay parameter at initial time s
%       ps_minus  - Left-tail decay parameter at initial time s
%
%   Outputs:
%       increments - An [N_sim x 1] vector of the simulated process increments
%
%   See also: POINTMASSCALCULATOR, SIMULATE_FROM_CDF

   if nargin < 14 || isempty(doPlot)
    doPlot = false;
    end

   if strcmp(activity, 'finite')
        % =========================================================================
        % STEP 1: CALCULATE THE JUMP PROBABILITY
        % =========================================================================
        % Retrieve Lambda (the total expected number of jumps)
        lambda = pointmasscalculator(pt_plus, pt_minus, ps_plus, ps_minus);
        
        % P(K > 0) = 1 - P(K = 0) = 1 - exp(-Lambda)
        p_jump = 1 - exp(-lambda);
    
        % =========================================================================
        % STEP 2: THE BERNOULLI COIN TOSS (MIXTURE MODEL SPLIT)
        % =========================================================================
        % Generate a boolean mask: 1 if the path jumped, 0 if it stayed flat.
        % Using 'rand < p' is heavily optimized in MATLAB for speed.
        jump_mask = rand(N_sim, 1) < p_jump;
        
        % Count exactly how many paths need complex CDF inversion
        N_jumps = sum(jump_mask);
    
        % =========================================================================
        % STEP 3: PREALLOCATE THE OUTPUT ARRAY
        % =========================================================================
        % Initialize all paths to 0 (the stochastic jump component)
        jump_component = zeros(N_sim, 1);
    
        % =========================================================================
        % STEP 4: CONDITIONAL JUMP SIMULATION
        % =========================================================================
        % Only trigger the heavy FFT machinery if at least one path jumped
        if N_jumps > 0
            % Call the simulation engine ONLY for the required number of paths
            [cdf_clean, x_grid_fine] = lewis_FFT_digital(@conditional_cf_MA_FA,M, dz, params,scale_factor, 1,'MA',doubleshift, doPlot);
            simulated_jumps = simulate_from_cdf(cdf_clean, x_grid_fine,spline,N_jumps);
            % Inject the simulated jumps back into the corresponding active paths
            jump_component(jump_mask) = simulated_jumps;
        end
        increments = jump_component + delta_mu;
        % NOTE: Paths where jump_mask == 0 naturally remain exactly 0 here.
   else 
        [cdf_clean, x_grid_fine] = lewis_FFT_digital(@cf_MA_IA,M, dz, params,scale_factor, 1,'MA',doubleshift, doPlot);
        increments = simulate_from_cdf(cdf_clean, x_grid_fine,spline,N_sim);
        
   end
 
end

