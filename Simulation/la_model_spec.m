function spec = la_model_spec(model, p)
%LA_MODEL_SPEC  Single point of divergence between the two "twin" Linear
%   Additive models with a Bachelier-type sqrt(T) scaling: Additive Bachelier
%   (AB) and Generalized Logistic (GL). Everything downstream (conditional CDF
%   via Lewis-FFT, tail adjustment, inverse-CDF simulation, forward-start MC) is
%   shared and model-agnostic: it only consumes the struct returned here.
%
%   spec = LA_MODEL_SPEC(MODEL, P)
%
%   INPUTS
%     model : 'AB' or 'GL' (case-insensitive)
%     p     : struct of CALIBRATED shape parameters
%               AB -> p.k   , p.eta
%               GL -> p.alpha, p.beta
%
%   OUTPUT (struct)
%     spec.model  : normalized model tag ('AB' / 'GL')
%     spec.cf     : marginal characteristic function handle of f_T,
%                   phi_T(u) = spec.cf(u, sigma, T), where the scale enters as
%                   the Bachelier dispersion sigma*sqrt(T). Accepts complex u.
%     spec.lamL   : left-tail exponential decay rate of the STANDARDIZED density
%                   (density ~ exp(+lamL*z) as z -> -inf). Sets the upper edge of
%                   the analyticity strip in the contour shift a (a < +lamL/sT).
%     spec.lamR   : right-tail exponential decay rate of the STANDARDIZED density
%                   (density ~ exp(-lamR*z) as z -> +inf). Sets the lower edge of
%                   the analyticity strip (a > -lamR/sT).
%     spec.varStd : variance of the STANDARDIZED increment (per unit
%                   (sigma*sqrt(T))^2). Used only to auto-size the spatial grid.
%     spec.I0     : normalization constant linking ATM vol to the model scale,
%                   sigma_t = sigma_ATM / I0  (Baviera-Massaria 2026, Eq. 14-15).
%
%   References:
%     - Baviera & Massaria (2026), Additive Bachelier (Biblio/3): AB CF, Eq. 4,
%       analyticity strip p+- = -+eta + sqrt(eta^2 + 1/k), I0 = sqrt(2pi) E[zeta_+].
%     - Carr & Torricelli, Additive logistic processes (Biblio/4): GL CF.

    switch upper(string(model))

        case "AB"
            k   = p.k;
            eta = p.eta;

            % Marginal CF of f_T at Bachelier scale sigma*sqrt(T): reuse the
            % single source of truth charateristic_function_AB (returns a handle).
            spec.cf = @(u, sigma, T) ab_cf(u, sigma, T, k, eta);

            % AB analyticity strip / tail rates (paper Eq. just below Eq. 4):
            %   right tail decay p+ , left tail decay p- , both > 0.
            p_plus  =  eta + sqrt(eta^2 + 1/k);
            p_minus = -eta + sqrt(eta^2 + 1/k);
            spec.lamR = p_plus;     % right tail  -> lower strip edge
            spec.lamL = p_minus;    % left  tail  -> upper strip edge

            % Var(f_T) = (1 + k*eta^2) * sigma_T^2 * T  => standardized variance:
            spec.varStd = 1 + k*eta^2;

            % I0 = sqrt(2*pi)*E[zeta_+] computed via Lewis-FFT at z=0 (I0.m).
            spec.I0 = I0(0, k, eta);

            spec.model = 'AB';

        case "GL"
            a = p.alpha;
            b = p.beta;

            % Marginal CF of f_T: phi_T(u) = cf_GL(u * sigma*sqrt(T)).
            spec.cf = @(u, sigma, T) cf_GL(u .* (sigma * sqrt(T)), a, b);

            % GL standardized density ~ exp(+a*z) (left) and exp(-b*z) (right),
            % so cf_GL is analytic for the shift a_shift in (-b/sT, +a/sT).
            spec.lamL = a;          % left  tail  -> upper strip edge
            spec.lamR = b;          % right tail  -> lower strip edge

            % Variance of the standardized Z-distribution = psi'(alpha)+psi'(beta).
            spec.varStd = psi(1, a) + psi(1, b);   % trigamma sum

            % I0 = sqrt(2*pi) * E[zeta_+] = sqrt(2*pi) * int_0^inf x f_std(x) dx
            % (Baviera-Massaria Eq. 14, same convention used in price_GL.m).
            integrand_mean = @(x) pdf_GL(a, b, x) .* x;
            spec.I0 = sqrt(2*pi) * quadgk(integrand_mean, 0, inf);

            spec.model = 'GL';

        otherwise
            error('la_model_spec:UnknownModel', ...
                  'model must be ''AB'' or ''GL'', got ''%s''.', string(model));
    end

    if ~isfinite(spec.I0) || spec.I0 == 0
        error('la_model_spec:I0Invalid', ...
              'I0 evaluated to a non-physical value (%g) for model %s.', ...
              spec.I0, spec.model);
    end
end

% -------------------------------------------------------------------------
function y = ab_cf(u, sigma, T, k, eta)
% Thin adapter so the AB CF matches the (u, sigma, T) signature: builds the
% handle once per call and evaluates it (charateristic_function_AB forces u(:)).
    phi = charateristic_function_AB(T, k, eta, sigma);
    y   = phi(u);
end
