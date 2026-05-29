function y = cf_MA_FA(u, params, scale_factor)
% CF_MA_FA Computes the characteristic function of the Finite Activity (FA)
% increment under the Minimal Additive (MA) model.
%
%   This function calculates the analytical Characteristic Function (CF) of the 
%   stochastic increment \Delta f = f_t - f_s. By the independent increments 
%   property of additive processes, this is exactly the ratio of the marginal 
%   characteristic functions at time t and s.
%
%   Formula:
%       \phi(u) = c * [ (ps_plus - iu)(ps_minus + iu) / (pt_plus - iu)(pt_minus + iu) ] * exp(i * \Delta\mu * u)
%       where c = (pt_plus * pt_minus) / (ps_plus * ps_minus)
%
%   Inputs:
%       u        : Frequency domain grid (can be a vector or matrix)
%       pt_plus  : Right tail parameter at target time t
%       pt_minus : Left tail parameter at target time t
%       ps_plus  : Right tail parameter at starting time s
%       ps_minus : Left tail parameter at starting time s
%       deltamu  : Deterministic drift of the increment (\mu_t - \mu_s)
%
%   Output:
%       y        : Evaluated characteristic function at u


    ps_plus = params(2)/scale_factor(1); ps_minus = params(1)/scale_factor(1);
    pt_plus = params(2)/scale_factor(2); pt_minus = params(1)/scale_factor(2);
    gamma_MA = (1 / params(1)) - (1 /params(2));
    deltamu = gamma_MA * (scale_factor(2) - scale_factor(1));

    % 1. Point mass probability 'c' (ratio of normalizing constants)
    c = (pt_plus * pt_minus) / (ps_plus * ps_minus);
    
    % 2. Complex rational fraction (numerator from s, denominator from t)
    numerator   = (ps_plus - 1i * u) .* (ps_minus + 1i * u);
    denominator = (pt_plus - 1i * u) .* (pt_minus + 1i * u);
    
    % 3. Deterministic drift phase shift
    phase_shift = exp(1i * deltamu .* u);
    
    % 4. Final CF assembly
    y = c .* (numerator ./ denominator) .* phase_shift;
    
end