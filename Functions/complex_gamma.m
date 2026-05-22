function y = complex_gamma(z)
% COMPLEX_GAMMA Computes the Gamma function for complex arguments
% using the Lanczos approximation (g=7, n=9) and Euler's reflection formula.
%
% Input:  z - Scalar, vector, or matrix of complex numbers
% Output: y - Evaluated Gamma function values with the same size as z

    % 1. Define Lanczos coefficients for g = 7, n = 9
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
     
    % Initialize the output array with the same dimensions as input
    y = zeros(size(z));
    
    % 2. Create logical masks to separate the two half-planes
    % Euler's reflection formula is required for Re(z) < 0.5
    reflect_mask = real(z) < 0.5;
    direct_mask  = ~reflect_mask;
    
    % 3. Branch 1: Direct Lanczos calculation (Re(z) >= 0.5)
    if any(direct_mask)
        z_dir = z(direct_mask);
        z_minus_1 = z_dir - 1;
        
        % Initialize the partial sum with the first coefficient
        x = p(1) * ones(size(z_dir));
        
        % Efficiently compute the rational function series
        for i = 2:length(p)
            x = x + p(i) ./ (z_minus_1 + (i-1));
        end
        
        % Compute the scaling factor
        t = z_minus_1 + g + 0.5;
        
        % Vectorized profile matching Stirling's asymptotic form
        y(direct_mask) = sqrt(2*pi) * (t.^(z_minus_1 + 0.5)) .* exp(-t) .* x;
    end
    
    % 4. Branch 2: Euler's Reflection Formula (Re(z) < 0.5)
    if any(reflect_mask)
        z_ref = z(reflect_mask);
        
        % Recursive call: since real(1 - z_ref) >= 0.5, the next step 
        % will naturally hit the direct branch, preventing infinite loops.
        y(reflect_mask) = pi ./ (sin(pi * z_ref) .* complex_gamma(1 - z_ref));
    end
    
    % 5. Numerical cleaning of infinitesimal imaginary parts
    % Removes floating-point noise on pure real numbers
    epsilon = 1e-12;
    pure_real_indices = abs(imag(y)) <= epsilon;
    y(pure_real_indices) = real(y(pure_real_indices));
end