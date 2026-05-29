function [moneyness_modified, c_mkt_calibration] = moneyness_generator(forward, strikes, calls, puts, sigma_ATM, yf, discount_factor, chi_min, chi_max)
% MONEYNESS_GENERATOR Computes modified moneyness and target call surface.
%
% Selects the most liquid OTM quotes per maturity within a band in the
% *normalized* moneyness chi := (K - F)/(sigma_ATM * sqrt(t)). OTM Puts
% are converted to synthetic Calls via Put-Call Parity. Strikes outside
% the band or lacking market data remain NaN.
%
% INPUTS:
%   forward         : (M x 1) forward prices per maturity
%   strikes         : (1 x N) market strike grid
%   calls           : (M x N) market call prices (NaNs allowed)
%   puts            : (M x N) market put prices  (NaNs allowed)
%   sigma_ATM       : (M x 1) Bachelier ATM implied vols
%   yf              : (M x 1) year fractions per maturity
%   discount_factor : (M x 1) discount factors per maturity
%   chi_min         : scalar lower bound on chi (e.g. -10), OTM put side
%   chi_max         : scalar upper bound on chi (e.g. +10), OTM call side
%
% OUTPUTS:
%   moneyness_modified : (M x N) chi values where market data is selected
%   c_mkt_calibration  : (M x N) target call prices (mid for OTM call, parity
%                        from put for OTM put)
%
% Rationale (intervento 2.7): a band in dollars is unfair across maturities
% because sigma * sqrt(t) grows with t; a band in chi covers a comparable
% number of standard deviations at every expiry, consistent with paper 3
% Sec. 3 (separability of implied vol in the moneyness degree).

    M = length(forward);
    N = length(strikes);

    moneyness_modified = NaN(M, N);
    c_mkt_calibration  = NaN(M, N);

    for i = 1:M
        curr_calls   = calls(i, :);
        curr_puts    = puts(i, :);
        curr_forward = forward(i);
        curr_yf      = yf(i);
        curr_df      = discount_factor(i);
        curr_sigma   = sigma_ATM(i);

        % Normalized moneyness chi = (K - F) / (sigma_ATM * sqrt(t)).
        % Same orientation as `strikes` (1 x N row).
        curr_chi = (strikes - curr_forward) / (curr_sigma * sqrt(curr_yf));

        % OTM Call leg: K >= F and chi <= chi_max.
        mask_otm_calls = (strikes >= curr_forward) & ...
                         (curr_chi <= chi_max)     & ...
                         ~isnan(curr_calls);

        % OTM Put leg: K < F and chi >= chi_min.
        mask_otm_puts  = (strikes <  curr_forward) & ...
                         (curr_chi >= chi_min)     & ...
                         ~isnan(curr_puts);

        % --- OTM calls: take mid as is ---
        c_mkt_calibration(i, mask_otm_calls)  = curr_calls(mask_otm_calls);
        moneyness_modified(i, mask_otm_calls) = curr_chi(mask_otm_calls);

        % --- OTM puts: synthetic call via parity C = P + DF*(F - K) ---
        derived_calls = curr_puts(mask_otm_puts) + ...
                        curr_df * (curr_forward - strikes(mask_otm_puts));
        c_mkt_calibration(i, mask_otm_puts)  = derived_calls;
        moneyness_modified(i, mask_otm_puts) = curr_chi(mask_otm_puts);
    end

end