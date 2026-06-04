function [moneyness_modified, c_mkt_calibration] = moneyness_generator(forward, strikes, calls, puts, sigma_ATM, yf, discount_factor, x_min, x_max)
% MONEYNESS_GENERATOR  OTM quote selection and target call surface ([3] Sec 4.1).
%   Selects OTM quotes per maturity in the dollar band x = K-F in [x_min, x_max]
%   ([-30,30]$). OTM puts become synthetic calls via parity C = P + DF*(F-K).
%   Returns the model coordinate chi = (K-F)/(sigma_ATM*sqrt(t)) ([3] Eq. 17);
%   strikes outside the band or without data stay NaN.
%
% INPUTS:
%   forward, sigma_ATM, yf, discount_factor : (M x 1) per maturity
%   strikes         : (1 x N) strike grid
%   calls, puts     : (M x N) market prices (NaNs allowed)
%   x_min, x_max    : dollar-moneyness band (OTM put / OTM call side)
% OUTPUTS:
%   moneyness_modified : (M x N) chi where data is selected
%   c_mkt_calibration  : (M x N) target call prices

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

        % Dollar moneyness x = K - F (selection band) and normalized
        % moneyness chi = x / (sigma_ATM * sqrt(t)) (model coordinate).
        % Both have the same orientation as `strikes` (1 x N row).
        curr_x   = strikes - curr_forward;
        curr_chi = curr_x / (curr_sigma * sqrt(curr_yf));

        % OTM Call leg: F <= K and x <= x_max.
        mask_otm_calls = (curr_x >= 0)      & ...
                         (curr_x <= x_max)  & ...
                         ~isnan(curr_calls);

        % OTM Put leg: K < F and x >= x_min.
        mask_otm_puts  = (curr_x <  0)      & ...
                         (curr_x >= x_min)  & ...
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