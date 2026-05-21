function sigma = sigmaATM(c_atm, B, yf, expiries)
% SIGMAATM  Bachelier ATM implied volatility from market call prices.
%
%   sigma = SIGMAATM(c_atm, B, yf) inverts the Bachelier call formula at
%   ATM (moneyness x = 0). At x = 0 the normalised Bachelier call reduces
%   to c_b(0,1) = phi(0) = 1/sqrt(2*pi), so the ATM call price is
%
%       C_ATM = B * sigma * sqrt(t) / sqrt(2*pi)
%
%   which gives the closed-form inversion (no iteration needed):
%
%       sigma^ATM = sqrt(2*pi / t) * C_ATM / B             [Eq. 15, 3]
%
%   Reference: Baviera & Massaria (2026), "The additive Bachelier model",
%   J. Comput. Appl. Math. 487, 117741, Eq. (15).
%
%   Inputs  (all column vectors, one entry per maturity)
%     c_atm    : ATM call prices   (output of callATM.m)
%     B        : discount factors  (output of bootstrap.m)
%     yf       : year fractions to each maturity (act/365)
%     expiries : (optional) datetime vector of expiry dates, used as x-axis labels
%
%   Output
%     sigma : Bachelier ATM implied volatilities (column vector)

c_atm = c_atm(:);
B     = B(:);
yf    = yf(:);

sigma = sqrt(2*pi ./ yf) .* (c_atm ./ B);

figure;
if nargin == 4
    plot(expiries, sigma, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
    xtickformat('MMM-yy');
    xlabel('Expiry');
else
    plot(yf, sigma, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
    xlabel('Time to maturity (yf)');
end
ylabel('\sigma^{ATM} ($/unit)');
title('Bachelier ATM Implied Volatility — Term Structure');
grid on;

end
