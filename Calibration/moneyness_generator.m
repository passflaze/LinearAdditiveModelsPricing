function [moneyness_modified, c_mkt_calibration] = moneyness_generator(forward, strikes, calls, puts, sigma_ATM, yf, discount_factor, x_min, x_max)
% MONEYNESS_GENERATOR Computes modified moneyness and target call surface.
%
% Selects the OTM quotes per maturity within a band in *dollar* moneyness
% x := K - F, as in Baviera & Massaria (2026), paper 3, Sec. 4.1 / Table 3
% (x in [-30 $, 30 $]). OTM Puts are converted to synthetic Calls via
% Put-Call Parity. Strikes outside the band or lacking market data remain
% NaN. The normalized moneyness chi := (K - F)/(sigma_ATM * sqrt(t)) is
% still returned as the model coordinate (it is what price_AB / price_MA /
% price_GL consume), but it no longer defines the selection band.
%
% INPUTS:
%   forward         : (M x 1) forward prices per maturity
%   strikes         : (1 x N) market strike grid
%   calls           : (M x N) market call prices (NaNs allowed)
%   puts            : (M x N) market put prices  (NaNs allowed)
%   sigma_ATM       : (M x 1) Bachelier ATM implied vols
%   yf              : (M x 1) year fractions per maturity
%   discount_factor : (M x 1) discount factors per maturity
%   x_min           : scalar lower bound on dollar moneyness (e.g. -30), OTM put side
%   x_max           : scalar upper bound on dollar moneyness (e.g. +30), OTM call side
%
% OUTPUTS:
%   moneyness_modified : (M x N) chi values where market data is selected
%   c_mkt_calibration  : (M x N) target call prices (mid for OTM call, parity
%                        from put for OTM put)
%
% Rationale: the dollar band [-30 $, 30 $] reproduces exactly the OTM set
% used in paper 3 (Table 3) and keeps this routine consistent with
% run_project2A (xMax = 30 $). chi remains the model coordinate via Eq. (17).

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