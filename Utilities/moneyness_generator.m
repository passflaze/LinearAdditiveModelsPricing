function [moneyness_modified, c_mkt_calibration_normed, norm_factor] = moneyness_generator(forward, strikes, calls, puts, sigma_ATM, yf, discount_factor, filter_type)
% MONEYNESS_GENERATOR Computes modified moneyness and normalized market prices.
%
% This function processes market option prices to construct the calibration 
% surface. It evaluates the modified moneyness grid and extracts the corresponding 
% call option prices. It incorporates a Put-Call Parity fallback for missing 
% call data and applies a structural filter to isolate specific moneyness regions.
%
% INPUTS:
%   forward         : Vector of forward prices for each maturity (length M)
%   strikes         : Vector of market strikes (length N)
%   calls           : Matrix (M x N) of market call prices (can contain NaNs)
%   puts            : Matrix (M x N) of market put prices (can contain NaNs)
%   sigma_ATM       : Vector of At-The-Money implied volatilities (length M)
%   yf              : Vector of year fractions / time to maturity (length M)
%   discount_factor : Vector of discount factors for each maturity (length M)
%   filter_type     : (Optional) String determining the moneyness filter. 
%                     'OTM' (Default) -> Only strikes > forward
%                     'ITM'           -> Only strikes < forward
%                     'ALL'           -> No filtering, uses all available strikes
%
% OUTPUTS:
%   moneyness_modified       : Matrix (M x N) of evaluated modified moneyness points
%   c_mkt_calibration_normed : Matrix (M x N) of target call prices scaled by the norm_factor
%   norm_factor              : Column vector (M x 1) used to scale the target prices

    % =========================================================================
    % STEP 1: INITIALIZATION & DEFAULT ARGUMENTS
    % =========================================================================
    % Set default filter to 'OTM' if the user does not provide the 8th argument
    if nargin < 8 || isempty(filter_type)
        filter_type = 'OTM';
    end

    M = length(forward);
    N = length(strikes);
    
    % Initialize empty matrices (unmatched slots will naturally remain NaN)
    moneyness_modified = NaN(M, N);
    c_mkt_calibration  = NaN(M, N);
    
    % =========================================================================
    % STEP 2: SURFACE POPULATION & FILTERING
    % =========================================================================
    for i = 1:M
        % 1. Extract row vectors for the current maturity row
        curr_calls = calls(i, :);
        curr_puts  = puts(i, :);
        
        % 2. Split filter into call-side and put-side masks.
        % For 'OTM': use OTM calls (K>F) AND OTM puts (K<F) via PCP -- both
        % sides of the smile contribute. Previously the same mask was applied
        % to puts as well, so the K<F branch of the surface was never seen by
        % the calibrator.
        switch upper(filter_type)
            case 'OTM'
                call_filter = strikes > forward(i);
                put_filter  = strikes < forward(i);
            case 'ITM'
                call_filter = strikes < forward(i);
                put_filter  = strikes > forward(i);
            case 'ALL'
                call_filter = true(1, N);
                put_filter  = true(1, N);
            otherwise
                error('MoneynessGenerator:InvalidFilter', 'filter_type must be ''OTM'', ''ITM'', or ''ALL''.');
        end

        % 3. Data-availability masks. For OTM/ITM the two filters are
        % disjoint by construction (K>F vs K<F); the ~mask_has_call guard
        % keeps 'ALL' safe by preferring direct call quotes over PCP.
        % Penny filter: drop quotes < 0.1 (unreliable / stale).
        penny_thr = 0.1;
        mask_has_call = ~isnan(curr_calls) & (curr_calls >= penny_thr) & call_filter;
        mask_put_only = ~isnan(curr_puts)  & (curr_puts  >= penny_thr) & put_filter & ~mask_has_call;
        
        % 4. Pre-compute the modified moneyness vector for the entire row
        curr_moneyness = (strikes - forward(i)) / (sigma_ATM(i) * sqrt(yf(i)));
        
        % --- CASE 1: Call option is directly available ---
        c_mkt_calibration(i, mask_has_call) = curr_calls(mask_has_call);
        moneyness_modified(i, mask_has_call) = curr_moneyness(mask_has_call);
        
        % --- CASE 2: Call is missing, but Put is available (Put-Call Parity) ---
        % Implied Call Price evaluation: C = P + DF * (F - K)
        derived_calls = curr_puts(mask_put_only) + discount_factor(i) * (forward(i) - strikes(mask_put_only));
        
        c_mkt_calibration(i, mask_put_only) = derived_calls;
        moneyness_modified(i, mask_put_only) = curr_moneyness(mask_put_only);
    end
    
    % =========================================================================
    % STEP 3: FINAL SCALING & NORMALIZATION
    % =========================================================================
    % Ensure structural parameters form a column vector using (:) to prevent 
    % dimensional broadcasting mismatch crashes during matrix/vector division.
    % Bachelier normalization (AB Eq. 5):  C/B_0 = sigma * sqrt(t) * c_b(...)
    % => norm_factor must be B_0 * sigma_ATM * sqrt(t), NOT B_0 * sigma_ATM * t.
    norm_factor = discount_factor(:) .* sigma_ATM(:) .* sqrt(yf(:));
    
    c_mkt_calibration_normed = c_mkt_calibration ./ norm_factor;

end