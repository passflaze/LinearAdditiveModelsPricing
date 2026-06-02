function [positions, hedgeGreeks] = build_hedge_AB(exoticGreeks, gC, gP, gF)
% BUILD_HEDGE_AB  Quantitative delta-gamma-vega hedging rule.
%
%   Finds the quantities (nC, nP, nF) of vanilla call, vanilla put and future
%   that neutralise the Delta, Gamma and Vega of +1 unit of the exotic.
%
%   The future carries delta = 1, gamma = 0, vega = 0, so the linear system
%   decouples:
%
%       [ gC.vega   gP.vega ] [ nC ]   [ exoticGreeks.vega  ]
%       [ gC.gamma  gP.gamma] [ nP ] = [ exoticGreeks.gamma ]
%
%       nF = exoticGreeks.delta - nC*gC.delta - nP*gP.delta
%
%   The book (long exotic, short hedge) is then approximately greek-flat.
%
% INPUTS:
%   exoticGreeks - (struct) .delta .gamma .vega  of the exotic (per 1 unit)
%   gC, gP, gF   - (struct) greeks of call, put, future (greeks_vanilla_AB)
%
% OUTPUTS:
%   positions   - (struct) .nC .nP .nF   quantities to HOLD to hedge +1 exotic
%   hedgeGreeks - (struct) residual greeks of the hedged book (should be ~0):
%                          .delta .gamma .vega

    % --- Solve the 2x2 (vega, gamma) block for [nC; nP] ------------------
    A = [gC.vega,  gP.vega; ...
         gC.gamma, gP.gamma];
    b = [exoticGreeks.vega; exoticGreeks.gamma];

    % Conditioning guard. A is singular when the call and put span the SAME
    % (vega, gamma) direction: this happens exactly when Kc == Kp, because by
    % put-call parity  P = C - (F-K)B  is linear in F and flat in sigma, so a
    % same-strike put has IDENTICAL gamma and vega as the call (only the delta
    % differs). Two strikes at different moneyness are therefore required.
    kappaA = cond(A);
    if kappaA > 1e8
        warning('build_hedge_AB:illConditioned', ...
            ['Hedge matrix nearly singular (cond = %.2e). Call and put carry ', ...
             'collinear (gamma, vega) — pick strikes at clearly different ', ...
             'moneyness (Kc ~= Kp). Falling back to pinv.'], kappaA);
        n = pinv(A) * b;
    else
        n = A \ b;
    end
    nC = n(1);
    nP = n(2);

    % --- Future closes the residual delta --------------------------------
    nF = exoticGreeks.delta - nC*gC.delta - nP*gP.delta;

    positions = struct('nC', nC, 'nP', nP, 'nF', nF);

    % --- Residual greeks of the hedged book (exotic - hedge) -------------
    hedgeGreeks = struct( ...
        'delta', exoticGreeks.delta - (nC*gC.delta + nP*gP.delta + nF*gF.delta), ...
        'gamma', exoticGreeks.gamma - (nC*gC.gamma + nP*gP.gamma + nF*gF.gamma), ...
        'vega',  exoticGreeks.vega  - (nC*gC.vega  + nP*gP.vega  + nF*gF.vega));
end
