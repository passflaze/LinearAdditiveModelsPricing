function c=callATM(C,P,K,forward)

% callATM: interpolate the call price at the forward (ATM) via cubic spline.
%   Applies the penny-option filter (price < 0.1) to both calls and puts
%   jointly before interpolating, following Azzone & Baviera (2021), Sec. 2.
% Inputs
%   C: call prices
%   P: put prices
%   K: strike prices
%   forward: forward price at which to interpolate the call price
% Output
%   c: interpolated call price at the forward   

mask = ~isnan(C) & ~isnan(P);
mask = mask & (C > 0.1) & (P > 0.1);

C = C(mask);
K = K(mask);


c= interp1(K,C,forward, "spline");


end
























