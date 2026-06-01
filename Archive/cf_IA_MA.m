function cf = cf_IA_MA(u, pt_plus, pt_minus,ps_plus,ps_minus,deltamu)
%CF_IA_MA Computes the Characteristic Function of the MA base process.
%
%   CF = CF_IA_MA(U, PT_PLUS, PT_MINUS) evaluates the exact 
%   characteristic function in the complex domain using only the tail 
%   decay parameters, dynamically computing the drift gamma_t.
%
%   Inputs:
%       u        - Frequency domain grid (vector or scalar)
%       pt_plus  - Right-tail decay parameter at time t
%       pt_minus - Left-tail decay parameter at time t
%
%   Outputs:
%       cf       - The evaluated Characteristic Function (complex array)

    % =========================================================================
    % 1. DYNAMIC PARAMETER EVALUATION
    % =========================================================================
    % Compute the deterministic drift gamma_t directly from the tail parameters
    gamma_t = (1 ./ pt_minus) - (1 ./ pt_plus);

    % =========================================================================
    % 2. POLYNOMIAL TERM EVALUATION
    % =========================================================================
    % Evaluates the argument of the logarithm:
    % 1 - i*u*(1/pt_plus - 1/pt_minus) + u^2 / (pt_plus * pt_minus)
    
    term_linear = 1i .* u .* ((1 ./ pt_plus) - (1 ./ pt_minus));
    term_quad   = (u.^2) ./ (pt_plus .* pt_minus);
    
    poly_term = 1 - term_linear + term_quad;
    
    % =========================================================================
    % 3. LOG-CHARACTERISTIC FUNCTION (LCF) & EXPONENTIATION
    % =========================================================================
    % Compute the LCF and transform it back to the CF domain
    
    lcf = -log(poly_term) + 1i .* gamma_t .* u;
    
    cf = exp(lcf);

end