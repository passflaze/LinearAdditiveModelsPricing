function y = complex_gamma(z)
% COMPLEX_GAMMA  Gamma function for complex arguments via Lanczos approximation
%   (g=7, n=9) and Euler's reflection formula for Re(z) < 0.5.
%
% INPUTS:
%   z - (scalar / vector / matrix) complex arguments
% OUTPUT:
%   y - Gamma(z), same size as z

    % Lanczos coefficients (g = 7, n = 9)
    g = 7;
    p = [0.99999999999980993, ...
         676.5203681218851, ...
         -1259.1392167224028, ...
         771.32342877765313, ...
         -176.61502916214059, ...
         12.507343278686905, ...
         -0.13857109526572012, ...
         9.9843695780195716e-6, ...
         1.5056327351493116e-7];

    y = zeros(size(z));

    reflect_mask = real(z) < 0.5;
    direct_mask  = ~reflect_mask;

    % Direct Lanczos path (Re(z) >= 0.5)
    if any(direct_mask)
        z_dir     = z(direct_mask);
        z_minus_1 = z_dir - 1;
        x = p(1) * ones(size(z_dir));
        for i = 2:length(p)
            x = x + p(i) ./ (z_minus_1 + (i-1));
        end
        t = z_minus_1 + g + 0.5;
        y(direct_mask) = sqrt(2*pi) * (t.^(z_minus_1 + 0.5)) .* exp(-t) .* x;
    end

    % Reflection path (Re(z) < 0.5) — Gamma(z) = pi / (sin(pi*z) * Gamma(1-z))
    if any(reflect_mask)
        z_ref  = z(reflect_mask);
        y_temp = pi ./ (sin(pi * z_ref) .* complex_gamma(1 - z_ref));

        % Gamma decays exponentially along the imaginary axis; NaN from Inf*0
        % at high imaginary frequencies is replaced by the correct limit 0.
        is_high_freq_nan = isnan(y_temp) & (abs(imag(z_ref)) > 10);
        y_temp(is_high_freq_nan) = 0;

        y(reflect_mask) = y_temp;
    end

    % Strip floating-point imaginary noise on pure-real inputs
    epsilon = 1e-12;
    pure_real = abs(imag(y)) <= epsilon;
    y(pure_real) = real(y(pure_real));

    if any(isnan(y(:)) & ~isnan(z(:)))
        warning('ComplexGamma:NaNOutput', ...
                'complex_gamma generated NaN values for valid inputs.');
    end
end