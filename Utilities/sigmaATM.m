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
%   Inputs  (all column vectors, one entry per maturity)
%     c_atm    : ATM call prices   (output of callATM.m)
%     B        : discount factors  (output of bootstrap.m)
%     yf       : year fractions to each maturity (act/365)
%     expiries : (optional) datetime vector of expiry dates, used as x-axis labels
%     doPlot   : (optional) boolean variables to govern the plot
%
%   Output
%     sigma : Bachelier ATM implied volatilities (column vector)


% --- Arguments Initialization ---
if nargin < 5 || isempty(doPlot)
    doPlot = false;
end
if nargin < 4 || isempty(expiries)
    expiries = [];
end
% -----------------------------------------

c_atm = c_atm(:);
B     = B(:);
yf    = yf(:);

sigma = sqrt(2*pi ./ yf) .* (c_atm ./ B);

% --- Plotting Block ---
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

