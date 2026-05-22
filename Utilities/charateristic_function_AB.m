function phi = charateristic_function_AB(t, k, eta, sigma_t)
% CHARACTERISTIC_FUNCTION_AB  Characteristic function of the Additive bacheleir model with alpha = 0.5

psi = @(z) (1./k) .* (1 - sqrt(1 + 2.*z.*k));

log_phi = @(u) psi(1i.*u.*(eta.*sigma_t.*sqrt(t)) + 0.5.*sigma_t.^2.*t.*u.^2) ...
              + 1i.*(eta.*sigma_t.*sqrt(t)).*u;

phi = @(u) exp(log_phi(u(:)));

end