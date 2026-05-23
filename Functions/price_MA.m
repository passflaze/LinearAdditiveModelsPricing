function CallPrice = price_MA(params, discount_factor, yearfraction, sigma_ATM, moneyness_modified)
%PRICE_MA Computes the Call Option price using the MA model.
%
% INPUTS:
%   params             : Struct containing model parameters (alpha, beta).
%                        Both alpha and beta should be strictly positive.
%   discount_factor    : Column vector (M x 1) of discount factors for each maturity.
%   yearfraction       : Column vector (M x 1) of time-to-maturities (in years).
%   sigma_ATM          : Column vector (M x 1) of ATM implied volatilities.
%   moneyness_modified : Matrix (M x K) of modified moneyness values, where M 
%                        is the number of maturities and K is the number of strikes.
%
% OUTPUT:
%   CallPrice          : Matrix (M x K) of calculated Call option prices.

    % 1. Extract parameters
    alpha = params(1);
    beta  = params(2);
    
    % 2. Calculate model constants
    gamma_MA = (1/alpha) - (1/beta);
    C     = 1 / ((1/beta) + (1/alpha)); % Equivalent to (1/beta + 1/alpha)^(-1)
    
    % 3. Calculate I0 based on the sign of gamma_MA (I0 is a scalar)
    I0 = I0_MA(gamma_MA,C,alpha,beta);
    
    % 4. Pre-compute the maturity-dependent multiplier
    % Since discount_factor, sigma_ATM, and yearfraction are all (M x 1),
    % this operation results in a column vector of size (M x 1).
    maturity_multiplier = discount_factor .* (sigma_ATM ./ I0) .* sqrt(yearfraction);
    
    % 5. Compute the exponential argument for the whole matrix
    % Resulting size is (M x K)
    X = I0 * moneyness_modified - gamma_MA;
    
    % 6. Vectorized conditional logic (replaces the scalar if-else)
    % Initialize the core exponential matrix with zeros (M x K)
    core_value = zeros(size(moneyness_modified));
    
    % Create logical masks for the two branches
    mask_less    = (X < 0);  % Equivalent to I0*moneyness_modified < gamma_MA
    mask_greater = ~mask_less; 
    
    % --- Branch 1: Condition is True ---
    core_value(mask_less) = (C / (alpha^2)) * exp(alpha * X(mask_less)) - moneyness_modified(mask_less)*I0;
    
    % --- Branch 2: Condition is False ---
    core_value(mask_greater) = (C / (beta^2)) * exp(-beta * X(mask_greater));
    
    % 7. Calculate final Call Price using implicit expansion
    CallPrice = maturity_multiplier .* core_value;

end