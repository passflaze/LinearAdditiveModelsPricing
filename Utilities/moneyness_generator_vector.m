function [moneyness_modified, c_mkt_calibration, yf_mapped, df_mapped, sigma_ATM_mapped] = moneyness_generator_vector(forward, strikes, calls, puts, sigma_ATM, yf, discount_factor, xmin, xmax)
% MONEYNESS_GENERATOR_VECTOR Computes 1D mapped vectors for calibration.
%
% This function flattens the options surface into 1D column vectors while 
% strictly selecting the most liquid market instruments within a defined 
% spatial boundary. It extracts OTM Calls (F <= K <= F + xmax) and OTM Puts 
% (F + xmin <= K < F), converting the latter into synthetic calls via Put-Call Parity. 
% Each valid option is then perfectly mapped to its corresponding year fraction, 
% discount factor, and ATM volatility.
%
% INPUTS:
%   forward         : Vector of forward prices for each maturity (length M)
%   strikes         : Vector of market strikes (length N)
%   calls           : Matrix (M x N) of market call prices (can contain NaNs)
%   puts            : Matrix (M x N) of market put prices (can contain NaNs)
%   sigma_ATM       : Vector of At-The-Money implied volatilities (length M)
%   yf              : Vector of year fractions (length M)
%   discount_factor : Vector of discount factors for each maturity (length M)
%   xmin            : Lower bound offset for Puts (e.g., -400)
%   xmax            : Upper bound offset for Calls (e.g., +400)
%
% OUTPUTS (All are K x 1 column vectors, where K is the total valid options):
%   moneyness_modified : Evaluated modified moneyness points (negative and positive)
%   c_mkt_calibration  : Target call prices (direct OTM and recovered from OTM Puts)
%   yf_mapped          : Year fraction mapped to each point
%   df_mapped          : Discount factor mapped to each point
%   sigma_ATM_mapped   : ATM volatility mapped to each point

    % =========================================================================
    % STEP 1: INITIALIZATION
    % =========================================================================
    [M, N] = size(calls);
    
    % Initialize empty column vectors for vertical accumulation
    moneyness_modified = [];
    c_mkt_calibration  = [];
    yf_mapped          = [];
    df_mapped          = [];
    sigma_ATM_mapped   = [];
    
    % =========================================================================
    % STEP 2: SURFACE POPULATION & SPATIAL FILTERING
    % =========================================================================
    for i = 1:M
        % 1. Extract row vectors and scalar values for the current maturity
        curr_calls   = calls(i, :);
        curr_puts    = puts(i, :);
        curr_forward = forward(i);
        curr_yf      = yf(i);
        curr_df      = discount_factor(i);
        curr_sigma   = sigma_ATM(i);
        
        % 2. Create logical masks combining Liquidity (OTM), Availability (~NaN), 
        %    and Spatial Boundaries (xmin, xmax).
        
        % OTM Calls: Forward <= Strike <= Forward + xmax
        mask_otm_calls = (strikes >= curr_forward) & ...
                         (strikes <= curr_forward + xmax) & ...
                         ~isnan(curr_calls);
        
        % OTM Puts: Forward + xmin <= Strike < Forward
        mask_otm_puts  = (strikes < curr_forward) & ...
                         (strikes >= curr_forward + xmin) & ...
                         ~isnan(curr_puts); 
        
        % 3. Pre-compute the modified moneyness vector for the entire row
        curr_moneyness = (strikes - curr_forward) / (curr_sigma * sqrt(curr_yf)); 
        
        % --- CASE 1: Assign filtered OTM Calls directly ---
        valid_calls_direct = curr_calls(mask_otm_calls);
        valid_money_direct = curr_moneyness(mask_otm_calls);
        
        % --- CASE 2: Recover filtered Calls from OTM Puts via Put-Call Parity ---
        % Synthetic Call = P + DF * (F - K)
        valid_calls_recovered = curr_puts(mask_otm_puts) + curr_df * (curr_forward - strikes(mask_otm_puts));
        valid_money_recovered = curr_moneyness(mask_otm_puts);
        
        % --- COMBINE DATA ---
        % Force everything to be column vectors using (:) and concatenate.
        combined_calls = [valid_calls_recovered(:); valid_calls_direct(:)];
        combined_money = [valid_money_recovered(:); valid_money_direct(:)];
        
        % Skip the rest of the loop if no valid options exist for this maturity
        num_valid = length(combined_calls);
        if num_valid == 0
            continue;
        end
        
        % --- ACCUMULATE INTO 1D MASTER VECTORS ---
        % Vertically concatenate the new data to the master arrays
        moneyness_modified = [moneyness_modified; combined_money];
        c_mkt_calibration  = [c_mkt_calibration; combined_calls];
        
        % Expand scalar mappings to match the precise number of valid points
        yf_mapped        = [yf_mapped; repmat(curr_yf, num_valid, 1)];
        df_mapped        = [df_mapped; repmat(curr_df, num_valid, 1)];
        sigma_ATM_mapped = [sigma_ATM_mapped; repmat(curr_sigma, num_valid, 1)];
    end
end