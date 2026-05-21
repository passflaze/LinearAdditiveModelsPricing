function c = callATM(C, P, K, forward, B)
% CALLATM  Interpolate the call price at the forward (ATM) via cubic spline.
%   Where call prices are missing (NaN), they are recovered from put prices
%   via put-call parity:  C = P + B*(F - K).
%   A penny-option filter (price < 0.1) is then applied to the puts used
%   for recovery, to avoid propagating unreliable quotes.
%   Points where both call and put are NaN are dropped entirely.
%
% Inputs
%   C       : call prices (row or column vector)
%   P       : put prices  (same size as C)
%   K       : strikes     (same size as C)
%   forward : ATM forward price (scalar)
%   B       : discount factor (scalar), required only when NaN calls exist
% Output
%   c : interpolated call price at the forward

C = C(:);  P = P(:);  K = K(:);

% recover missing calls from puts via put-call parity (penny filter on puts)
maskNanC = isnan(C);
if any(maskNanC)
    validPut = maskNanC & ~isnan(P) & (P > 0.1);
    C(validPut) = P(validPut) + B * (forward - K(validPut));
end

% drop strikes where call is still NaN after recovery
mask = ~isnan(C);
K = K(mask);
C = C(mask);

c = interp1(K, C, forward, 'spline');


end
























