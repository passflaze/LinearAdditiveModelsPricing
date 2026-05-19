function [discount_factor, forward] = bootstrap(P, C, K)
% BOOTSTRAP  Synthetic discount factor and forward price from option prices.
%
%   [discount_factor, forward] = BOOTSTRAP(P, C, K) implements the
%   synthetic-forward technique of Azzone & Baviera (2022), "Synthetic
%   forwards and cost of funding in the equity derivative market", for a
%   single maturity T. Starting from put-call parity
%
%       C_i - P_i = D(t,T) * ( F(t,T) - K_i ),
%
%   the function regresses G_i := C_i - P_i on K_i via OLS:
%
%       G = beta0 + beta1 * K,    with   beta0 = D*F,  beta1 = -D,
%
%   so that the discount factor is D = -beta1 and the forward price is
%   F = beta0 / D. Vectors P, C, K must contain mid-quotes of European
%   puts, calls and the corresponding strikes for the same maturity.
%
%   Inputs are assumed to be mid-quotes. Only the penny-option filter
%   (price > 0.1) is applied; the bid-ask spread filter of the paper
%   is not, as bid/ask data are not available.
%
%   Inputs
%     P : put mid-prices  (vector)
%     C : call mid-prices (vector, same length as P)
%     K : strikes         (vector, same length as P)
%
%   Outputs
%     discount_factor : scalar D(t,T)
%     forward         : scalar F(t,T)

if length(P) ~= length(C) || length(P) ~= length(K)
    error("P, C and K must have the same length");
end

% force column vectors so the design matrix is well-shaped regardless
% of how the caller passed the inputs
P = P(:);
C = C(:);
K = K(:);

% penny-option filter only: both legs must be quoted and above 0.1
% (bid-ask spread filter from the paper not applicable: only mid-quotes available)
mask = ~isnan(C) & ~isnan(P);
mask = mask & (C > 0.1) & (P > 0.1);

C = C(mask);
P = P(mask);
K = K(mask);

% OLS on put-call parity: G = -D*K + D*F  =>  x(1) = -D, x(2) = D*F
G = C - P;
A = [K, ones(length(K), 1)];
x = A \ G;

discount_factor = -x(1);
forward         =  x(2) / discount_factor;

end