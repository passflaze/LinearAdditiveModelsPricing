function cf = cf_MA_IA(u, params, scale_factors)
%CF_MA_IA Characteristic function of the MA marginal (base process at time t).
%
%   CF = CF_MA_IA(U, PARAMS, SCALE_FACTORS) evaluates the exact marginal CF of
%   the Minimal Additive process. The tail decays are derived from the model
%   parameters and the (largest) scale factor, pt_plus = beta/scale,
%   pt_minus = alpha/scale, and the centering drift gamma_t is computed so that
%   the marginal has zero mean (martingale).
%
%   Inputs:
%       u             - Fourier argument (vector or scalar)
%       params        - [alpha; beta] MA shape parameters
%       scale_factors - scale factor(s); the largest entry is used,
%                       scale = max(scale_factors) = sigma_t*sqrt(t)
%
%   Outputs:
%       cf       - The evaluated characteristic function (complex array)

    scale_factor = max(scale_factors);
    pt_plus = params(2)/scale_factor; pt_minus = params(1)/scale_factor;
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