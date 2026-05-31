function sigma = bachelier_iv(C, B, fwdMinusK, t)
% BACHELIER_IV  Invert (discounted) Bachelier call prices to the implied vol.
%
%   sigma = BACHELIER_IV(C, B, fwdMinusK, t) solves, via Newton, for the
%   Bachelier (normal) implied volatility sigma in
%
%       C = B * [ (F-K) * Phi(d) + v * phi(d) ],   d = (F-K)/v,  v = sigma*sqrt(t)
%
%   where Phi/phi are the standard normal CDF/PDF. The map is monotone and
%   smooth in v (vega dC/dv = B*phi(d) > 0), so Newton converges fast.
%   Vectorized over C and fwdMinusK; non-converged/degenerate inversions are
%   returned as NaN.
%
% INPUTS
%   C         : (vector) discounted call price(s)
%   B         : discount factor(s); scalar (broadcast) OR per-element (size of C)
%   fwdMinusK : (vector) F - K  (dollar moneyness), same size as C
%   t         : year fraction(s); scalar (broadcast) OR per-element (size of C)
%
%   The per-element form for B and t lets the WHOLE (M x N) implied-vol surface
%   be inverted in a single vectorized call (flatten C, fwdMinusK and the
%   per-maturity B, t), instead of looping over maturities.
%
% OUTPUT
%   sigma     : (column vector) Bachelier implied vol(s) (NaN where it failed)

    C = C(:); fwdMinusK = fwdMinusK(:);
    n = numel(C);

    % Broadcast scalar discount/maturity to per-element column vectors.
    B = B(:); t = t(:);
    if isscalar(B), B = repmat(B, n, 1); end
    if isscalar(t), t = repmat(t, n, 1); end

    sigma = nan(n, 1);

    % ATM-style starting guess: at F=K, C = B*v/sqrt(2*pi) => v0 = C/B*sqrt(2pi).
    v = max(C ./ B * sqrt(2*pi), 1e-8);

    for it = 1:100
        d   = fwdMinusK ./ v;
        Phi = 0.5 * erfc(-d / sqrt(2));
        phi = exp(-0.5 * d.^2) / sqrt(2*pi);
        f   = B .* (fwdMinusK .* Phi + v .* phi) - C;  % residual
        vega = B .* phi;                               % dC/dv
        vega(vega < 1e-14) = 1e-14;                    % guard flat regions
        step = f ./ vega;
        v = v - step;
        v = max(v, 1e-10);                             % keep positive
        if max(abs(step)) < 1e-12, break; end
    end

    % Discard non-converged / degenerate inversions.
    d   = fwdMinusK ./ v;
    Phi = 0.5 * erfc(-d / sqrt(2));
    phi = exp(-0.5 * d.^2) / sqrt(2*pi);
    resid = B .* (fwdMinusK .* Phi + v .* phi) - C;
    ok = abs(resid) < 1e-6 * max(C, 1e-8);

    sigma(ok) = v(ok) ./ sqrt(t(ok));
end
