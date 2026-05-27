function [k,eta,sigma_t,rMSE] = calibrateAB(chi, cMkt, normFact, sigmaATM)
% CALIBRATEAB  Calibrate the Additive Bachelier model parameters (k, eta).
%
% Follows Baviera & Massaria (2026), Eq. (20): the loss is the L^2 distance
% between MARKET CALL PRICES IN $ and MODEL CALL PRICES IN $
%
%     d(eta, k) = sum_i [ C^mkt_i  -  B_i * sigma^ATM_i * sqrt(t_i) * G(chi_i; eta, k) ]^2
%
% where G(chi; eta, k) is the normalized model price (Eq. 19), computed by
% call_AB_FFT (which is built so that call_AB_FFT(chi, k, eta) = G(chi; eta, k)).
%
% Inputs (column vectors, one entry per quoted option pooled across maturities)
%   chi      : normalized moneyness chi = (K - F) / (sigma_ATM * sqrt(t))   (Eq. 17)
%   cMkt     : market call price in $ (OTM calls direct, OTM puts via put-call parity)
%   normFact : B_0 * sigma_ATM * sqrt(t) for each point  (same maturity -> same value)
%   sigmaATM : ATM Bachelier vols, one per maturity (used only to recover sigma_t at the end)
%
% Outputs
%   k, eta  : calibrated AB parameters
%   sigma_t : term structure of the model scale, sigma_t = sigma_ATM / I_0   (Eq. 15)
%   rMSE    : root mean squared error in $ on the calibration set

chi      = chi(:);
cMkt     = cMkt(:);
normFact = normFact(:);

% Model price in $ via Eq. (20):  C^mod = (B * sigma_ATM * sqrt(t)) * G(chi; eta, k)
% I_0 is computed ONCE per (k, eta) trial and reused inside call_AB_FFT
% (otherwise call_AB_FFT recomputes it via its own 2^14 FFT, doubling work).
cMod = @(k_p, eta_p) normFact .* call_AB_FFT(chi, k_p, eta_p, I0(0, k_p, eta_p));

% L^2 distance on full $ prices. Wrap with safeObj (local function) to penalize
% non-finite outputs so fmincon steps back into the feasible region instead of crashing.
objective = @(p) safeObj(p, cMod, cMkt);

% Bounds chosen from paper Fig.2-3 (Covid period: eta ~ [-0.3, 0.3], k ~ [0.4, 1.2]),
% with generous margin and far from the strip boundary where the Lewis contour blows up.
initial_guess = [1.0, 0.2];
lb = [1e-3, -1.5];
ub = [ 5.0,  1.5];

options = optimoptions('fmincon', 'Display', 'iter', ...
                       'Algorithm', 'interior-point', ...
                       'OptimalityTolerance', 1e-8, ...
                       'StepTolerance',       1e-10, ...
                       'MaxFunctionEvaluations', 5000);

[param_opt, SSE] = fmincon(objective, initial_guess, [], [], [], [], ...
                           lb, ub, [], options);

k   = param_opt(1);
eta = param_opt(2);

% I_0 = sqrt(2*pi) * E[c_b(eta*(G-1), sqrt(G))]   (Eq. 14)
% Computed via Lewis FFT in I0.m (sigma_t = 1 inside the CF, so the FFT at z=0 returns I_0).
I_0     = I0(0, k, eta);
sigma_t = sigmaATM(:) / I_0;          % Eq. (15):  sigma_ATM = sigma_t * I_0

rMSE = sqrt(SSE / numel(cMkt));       % root mean squared error in $

end

% --- local function ---------------------------------------------------------
function val = safeObj(p, cMod, cMkt)
% Robust SSE wrapper: any throw or non-finite output is mapped to a huge value,
% so fmincon backtracks instead of crashing on interp1 with all-NaN data.
    try
        r = cMod(p(1), p(2)) - cMkt;
    catch
        val = 1e20;
        return
    end
    if any(~isfinite(r))
        val = 1e20;
    else
        val = sum(r.^2);
    end
end






















