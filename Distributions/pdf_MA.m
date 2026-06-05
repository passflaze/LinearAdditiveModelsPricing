function pdf_val = pdf_MA(params, scale_factor, x)
% PDF_MA  Probability density function of the MA stochastic increment.
%   Asymmetric Laplace (double-exponential) centred at deltamu = gamma_MA * scale.
%
% INPUTS:
%   params       - [alpha; beta] left/right tail decay parameters (> 0)
%   scale_factor - integrated volatility sigma * sqrt(t); max() is taken
%   x            - (vector) evaluation points
% OUTPUT:
%   pdf_val      - (vector) density values at x

    alpha        = params(1);
    beta         = params(2);
    scale_factor = max(scale_factor);
    pt_plus      = beta  / scale_factor;
    pt_minus     = alpha / scale_factor;

    A       = (pt_plus * pt_minus) / (pt_plus + pt_minus);
    gamma_MA = (1 / alpha) - (1 / beta);
    deltamu  = gamma_MA * scale_factor;

    pdf_val = zeros(size(x));

    idx_left  = x < deltamu;
    idx_right = x >= deltamu;

    if any(idx_left)
        pdf_val(idx_left)  = A * exp( pt_minus * (x(idx_left)  - deltamu));
    end
    if any(idx_right)
        pdf_val(idx_right) = A * exp(-pt_plus  * (x(idx_right) - deltamu));
    end

end