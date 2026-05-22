function callPrice = europeanCallLA(F, B, K, phi, a)
%  !!!TO BE DONE WITH FFT!!!
% EUROPEANCALLLA  European Call price for Linear Additive processes.
%.  
%   callPrice = EUROPEANCALLLA(F, B, K, phi)
%   callPrice = EUROPEANCALLLA(F, B, K, phi, a)
%
%   Implements the European Call pricing formula for a Linear Additive
%   process, i.e. the "Call option Bachelier-Lewis formula" of
%   Proposition 2.2, eq. (8), in Baviera & Massaria (2026), "The Additive
%   Bachelier model with an application to the oil option market in the
%   Covid period", JCAM 487, 117741 (reference [3] of the project). This
%   covers point 2.a.
%
%   In a Linear Additive model the forward at maturity t is
%
%       F_t = F0 + f_t ,
%
%   where {f_t} is an additive, martingale driving process (E[f_t] = 0)
%   with characteristic function phi_t(u) = E0[ exp(i u f_t) ], analytic
%   in the horizontal strip Im(u) in (-p_t^+, p_t^-), p_t^+- > 0. The
%   Call with strike K and moneyness x := K - F0 has discounted price
%   C = B0 * E0[ F_t - K ]^+ = B0 * E0[ f_t - x ]^+, which by eq. (8) of
%   [3] equals, for any contour offset a in (-p_t^+, p_t^-),
%
%       C = B0 * ( R_a + (e^{x a} / 2pi) *
%               Int_{-Inf}^{Inf} [ phi_t(xi+i a) e^{-i xi x} / (i xi - a)^2
%                                  + 1{a=0} / xi^2 ] d xi ) ,
%
%   with the residue term
%
%       R_a = 0      if a in (-p_t^+, 0)
%       R_a = -x/2   if a = 0
%       R_a = -x     if a in (0, p_t^-) .
%
%   Since the price is real and the integrand at -xi is the conjugate of
%   the integrand at xi, the integral is computed as 2*Int_0^Inf of the
%   real part: this also cancels the odd 1/xi singularity that appears in
%   the a = 0 case, leaving a numerically smooth integrand.
%
%   Inputs
%     F   : scalar forward F0 = F(t0,t) for the maturity considered
%     B   : scalar discount factor B0 = B(t0,t) for the same maturity
%     K   : strikes (vector) -- e.g. the Fwd-ITM grid of point 2)
%     phi : function handle, phi(u) = characteristic function of the
%           driving process f_t. It MUST accept a complex argument, since
%           the formula evaluates phi_t on the line Im(u) = a.
%     a   : (optional) contour offset, must lie in the analyticity strip
%           (-p_t^+, p_t^-). Default a = 0, which always lies in the
%           strip and needs no knowledge of the tail coefficients; a ~= 0
%           gives a fully analytic (singularity-free) integrand and is
%           the choice generally used in practice (see [3], Sec. 2.1).
%
%   Output
%     callPrice : Call prices, same size as K (column vector).

if nargin < 5 || isempty(a)
    a = 0;                       % default contour: a = 0 (always in strip)
end

% column vector so the formula is shape-stable regardless of the caller
K = K(:);
x = K - F;                       % moneyness x := K - F0, one per strike

% residue term R_a (eq. 8): depends on the sign of the contour offset a
if a == 0
    Ra = -x / 2;
elseif a > 0
    Ra = -x;                     % a in (0, p_t^-)
else
    Ra = zeros(size(x));         % a in (-p_t^+, 0)
end

% I = Int_0^Inf Re[ integrand ] d xi, one value per strike.
% xi is scalar inside integral(...), x is the vector of moneynesses, so
% the integrand returns one row per evaluation point (ArrayValued).
integrand = @(xi) integrandLA(xi, x, phi, a);

I = integral(integrand, 0, Inf, ...
             'ArrayValued', true, 'RelTol', 1e-10, 'AbsTol', 1e-12);

% C = B0 * ( R_a + e^{x a} / (2 pi) * Int_{-Inf}^{Inf} = R_a + e^{x a}/pi * Int_0^Inf Re )
callPrice = B .* ( Ra + exp(x .* a) .* I(:) / pi );

end

% -------------------------------------------------------------------------

function val = integrandLA(xi, x, phi, a)
% Real part of the eq. (8) integrand at frequency xi, for every moneyness:
%
%   Re[ phi_t(xi + i a) e^{-i xi x} / (i xi - a)^2 + 1{a=0} / xi^2 ] .
%
% For a ~= 0 the integrand is analytic everywhere, xi = 0 included. For
% a = 0 the two terms combine to (1 - phi_t(xi) e^{-i xi x}) / xi^2, whose
% real part has a removable singularity at xi = 0 (phi is the cf of a
% centered rv): there the limiting value Re = (x^2 + Var(f_t)) / 2 is
% recovered with a small finite-difference probe.

if a == 0 && xi == 0
    dxi  = 1e-6;                                   % one-sided probe
    base = real( (1 - phi(dxi) .* exp(-1i*dxi*x)) ./ dxi^2 );
    val  = base;
    return;
end

core = phi(xi + 1i*a) .* exp(-1i*xi*x) ./ (1i*xi - a).^2;

if a == 0
    core = core + 1 ./ xi.^2;                      % 1{a=0}/xi^2 correction
end

val = real(core);

end
