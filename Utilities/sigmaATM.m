function sigma = sigmaATM(c_atm, B, yf, expiries, doPlot)
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
% INPUTS:
%   c_atm    - (column vector) ATM call prices, one per maturity
%   B        - (column vector) discount factors
%   yf       - (column vector) year fractions (act/365)
%   expiries - (column vector, optional) datetime expiry dates; used as x-axis labels
%   doPlot   - (logical, optional, default false) draw the term-structure figure
%
% OUTPUT:
%   sigma    - (column vector) Bachelier ATM implied volatilities

if nargin < 5 || isempty(doPlot)
    doPlot = false;
end
if nargin < 4 || isempty(expiries)
    expiries = [];
end

c_atm = c_atm(:);
B     = B(:);
yf    = yf(:);

sigma = sqrt(2*pi ./ yf) .* (c_atm ./ B);

    if doPlot
        figure('Name', 'Bachelier ATM Implied Volatility', 'Color', 'white');
        
        if ~isempty(expiries)
            plot(expiries, sigma, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
            xtickformat('MMM-yy');
            xlabel('Expiry', 'Interpreter', 'latex');
        else
            plot(yf, sigma, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
            xlabel('Time to maturity (yf)', 'Interpreter', 'latex');
        end
        
        ylabel('$\sigma^{ATM}$ (\$/unit)', 'Interpreter', 'latex');
        title('Bachelier ATM Implied Volatility — Term Structure', 'Interpreter', 'latex');
        grid on;
    end
end

