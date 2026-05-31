function pdf_val = pdf_MA(params, scale_factor, x)
% PDF_MA Computes the Probability Density Function (PDF) of the Minimal Additive
% stochastic increment at given point(s) x.
%
% Inputs:
%   alpha  - Left tail decay parameter (> 0)
%   beta   - Right tail decay parameter (> 0)
%   scale_factor - Integrated volatility (sigma * sqrt(t))
%   x      - Evaluation point(s) for the density (scalar or array)
%
% Output:
%   pdf_val - Evaluated density at x

    % 1. Calculate the scaled tail decay parameters
    alpha = params(1);
    beta = params(2);
    scale_factor = max(scale_factor);
    pt_plus  = beta / scale_factor;
    pt_minus = alpha / scale_factor;

    % 2. Calculate normalization constant A
    % A = (1/pt_plus + 1/pt_minus)^(-1)
    A = (pt_plus * pt_minus) / (pt_plus + pt_minus);

    % 3. Calculate the shift (deltamu)
    gamma_MA = (1 / alpha) - (1 / beta);
    deltamu = gamma_MA * scale_factor;

    % 4. Initialize output array to match the size of x
    pdf_val = zeros(size(x));

    % 5. Create logical masks for the piecewise evaluation
    idx_left  = x < deltamu;
    idx_right = x >= deltamu;

    % 6. Evaluate the PDF on the Left Branch
    % Formula: A * exp(pt_minus * (x - deltamu))
    if any(idx_left)
        pdf_val(idx_left) = A * exp(pt_minus * (x(idx_left) - deltamu));
    end

    % 7. Evaluate the PDF on the Right Branch
    % Formula: A * exp(-pt_plus * (x - deltamu))
    if any(idx_right)
        pdf_val(idx_right) = A * exp(-pt_plus * (x(idx_right) - deltamu));
    end

end