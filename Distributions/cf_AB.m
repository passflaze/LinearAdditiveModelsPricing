function phi = cf_AB(u, params, scale_factor)
% CF_AB  Computes the characteristic function of the log-price process
%        under the Additive Bachelier (AB) model with stability index alpha = 0.5,
%        corresponding to a tempered stable (inverse Gaussian) distribution.
%
%        The log-CF is given by:
%
%          ln phi_s(u) = psi( i*u*eta*scale_factor + 0.5*(u*scale_factor)^2; kappa, 0.5 )
%                        + i*u*eta*scale_factor
%
%        where scale_factor = sigma_s * sqrt(s) absorbs the volatility and
%        time scaling, and psi is the Laplace exponent of a tempered stable
%        random variable evaluated at alpha = 0.5:
%
%          psi(z; kappa, 0.5) = (1/kappa) * (1 - sqrt(1 + 2*z*kappa))
%
% INPUTS:
%   u            - (vector) Fourier argument (real-valued evaluation points)
%   params       - (2x1 vector) Model parameters, project-wide [k; eta]:
%                    params(1) = k   : mean-reversion / tempering parameter (k_s)
%                    params(2) = eta : drift parameter (eta_s)
%   scale_factor - (scalar) Composite scaling factor sigma_s * sqrt(s),
%                    where sigma_s is the volatility and s is the time horizon
%
% OUTPUT:
%   phi          - (vector, same size as u) Values of the characteristic
%                    function evaluated at u

    % Model parameters (project-wide convention [k; eta])
    k   = params(1);
    eta = params(2);

    % Laplace exponent (alpha = 0.5)
    psi = @(z) (1./k) .* (1 - sqrt(1 + 2.*z.*k));

    z_arg = 1i .* u .* (eta .* scale_factor) + 0.5 .* (scale_factor .* u).^2;

    log_phi = psi(z_arg) + 1i .* (eta .* scale_factor) .* u;
    
    phi = exp(log_phi);

end