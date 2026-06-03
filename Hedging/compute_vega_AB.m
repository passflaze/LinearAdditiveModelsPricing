function vega = compute_vega_AB(type, mkt, params_hedge, mc, bump, maturity_index, K)
% COMPUTE_VEGA_AB Computes the recalibrated Vega for exotic or vanilla options under the AB model.
%
%   This function calculates the sensitivity of the option price to a parallel
%   shift in the ATM market volatility surface (Vega). It uses a central finite
%   difference scheme. Crucially, the model is recalibrated for both the "up"
%   and "down" shifts to capture the full recalibration effect (Shadow Vega).
%
% INPUTS:
%   type           - (string) Option type identifier (e.g., 'CoC', 'PoP', 'Chooser', 'vanilla')
%   mkt            - (struct) Market data containing sigma_ATM, yf, discount_factor, etc.
%   params_hedge   - (struct) parameters for hedging/pricing. For vanillas, must contain .Kc and .forward
%   mc             - (struct) Monte Carlo/FFT numerical parameters (M, dz)
%   bump           - (scalar or vector) The shock size applied to the volatility surface
%   maturity_index - (scalar, optional) maturity leg for the VANILLA branch
%                    (default 4 = T2). Ignored by the exotic branch, which
%                    always uses the [iT1, iT2] = [2, 4] legs.
%   K              - (scalar, optional) strike for the VANILLA branch
%                    (default params_hedge.Kc). Ignored by the exotic branch.
%
% OUTPUTS:
%   vega         - (scalar) The estimated Vega of the option.

    if nargin < 6 || isempty(maturity_index)
        maturity_index = 4;   % default: T2 leg
    end
    if nargin < 7 || isempty(K)
        K = params_hedge.Kc;  % default strike (back-compat)
    end

    % =========================================================================
    % STEP 1: APPLY BUMP TO VOLATILITY SURFACE (PARALLEL SHIFT)
    % =========================================================================
    sigma_ATM_bumped_up   = mkt.sigma_ATM + bump;
    sigma_ATM_bumped_down = mkt.sigma_ATM - bump;
    
    % Setup optimization options for fmincon (suppress output to keep console clean)
    options_AB = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp');
    
    % Initial guess and bounds for AB calibration
    x0_AB = [1.0, 0.2];
    lb_AB = [1e-3, -1.5];
    ub_AB = [5.0,   1.5];
    
    % =========================================================================
    % STEP 2: RECALIBRATION
    % =========================================================================
    % 2A. Calibrate to UP bumped surface
    obj_fun_AB_up = @(x) objective_function_AB(x, mkt.discount_factor, mkt.yf, ...
                         sigma_ATM_bumped_up, mkt.moneyness_modified, mkt.c_mkt_calibration);
                   
    [params_AB_vega_up, ~, ~] = fmincon(obj_fun_AB_up, x0_AB, [], [], [], [], lb_AB, ub_AB, [], options_AB);
    
    % 2B. Calibrate to DOWN bumped surface
    obj_fun_AB_down = @(x) objective_function_AB(x, mkt.discount_factor, mkt.yf, ...
                           sigma_ATM_bumped_down, mkt.moneyness_modified, mkt.c_mkt_calibration);
                   
    [params_AB_vega_down, ~, ~] = fmincon(obj_fun_AB_down, x0_AB, [], [], [], [], lb_AB, ub_AB, [], options_AB);
    
    % =========================================================================
    % STEP 3: REPRICING
    % =========================================================================
    % Recompute scale factors using the new parameters and bumped ATM vols.
    scale_factor_up = [ (sigma_ATM_bumped_up(2) / I0_AB(0, params_AB_vega_up)) * sqrt(mkt.yf(2)), ...
                        (sigma_ATM_bumped_up(4) / I0_AB(0, params_AB_vega_up)) * sqrt(mkt.yf(4)) ];
                        
    scale_factor_down = [ (sigma_ATM_bumped_down(2) / I0_AB(0, params_AB_vega_down)) * sqrt(mkt.yf(2)), ...
                          (sigma_ATM_bumped_down(4) / I0_AB(0, params_AB_vega_down)) * sqrt(mkt.yf(4)) ];
    
    if strcmp(type, 'vanilla')
        % VANILLA PRICING
        % Strike K is passed explicitly (fixes the old Kc-for-puts bug); the
        % forward is read at the vanilla's own maturity leg, consistent with
        % greeks_vanilla_AB (which uses mkt.forward(maturity_index)).
        F = mkt.forward(maturity_index);

        % Scale factor recomputed at the vanilla's OWN maturity leg, using the
        % recalibrated parameters and the bumped ATM vols (no hardcoded leg).
        mi              = maturity_index;
        sf_up_vanilla   = (sigma_ATM_bumped_up(mi)   / I0_AB(0, params_AB_vega_up))   * sqrt(mkt.yf(mi));
        sf_down_vanilla = (sigma_ATM_bumped_down(mi) / I0_AB(0, params_AB_vega_down)) * sqrt(mkt.yf(mi));

        % Discount to t0 (same PV measure as greeks_vanilla_AB): lewis_FFT_call
        % returns the undiscounted forward call value, hence the B factor.
        B = mkt.discount_factor(mi);
        price_up   = B * lewis_FFT_call(@cf_AB, mc.M, mc.dz, params_AB_vega_up,   sf_up_vanilla,   K - F, 1, 'AB');
        price_down = B * lewis_FFT_call(@cf_AB, mc.M, mc.dz, params_AB_vega_down, sf_down_vanilla, K - F, 1, 'AB');
    else
        % EXOTIC PRICING
        price_up   = price_exotic_AB(type, params_AB_vega_up, scale_factor_up, mkt, params_hedge, mc);
        price_down = price_exotic_AB(type, params_AB_vega_down, scale_factor_down, mkt, params_hedge, mc);
    end
    
    % =========================================================================
    % STEP 4: CENTRAL DIFFERENCES AND VEGA COMPUTATION
    % =========================================================================
    % Calculate Vega using the central difference quotient. 
    % If bump is a vector, use bump(1) to ensure the denominator is a scalar.
    vega = (price_up - price_down) / (2 * bump(1));
end