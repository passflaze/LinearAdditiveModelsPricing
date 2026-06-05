function lambda = pointmasscalculator(pt_plus, pt_minus, ps_plus, ps_minus)
%POINTMASSCALCULATOR Computes the total expected number of jumps (Lambda).
%   LAMBDA = POINTMASSCALCULATOR(PT_PLUS, PT_MINUS, PS_PLUS, PS_MINUS) 
%   calculates the integrated intensity of the Lévy jump measure, Lambda, 
%   for a Minimal Additive (MA) process increment over the time interval [s, t].
%
%   In the MA framework, the probability of observing zero jumps (the size 
%   of the probability atom/point mass at the martingale-corrected drift) 
%   is given by c = exp(-lambda) = (pt_plus * pt_minus) / (ps_plus * ps_minus).
%   Therefore, lambda represents the negative logarithm of the atom size.
%
%   Inputs:
%       pt_plus  - Right-tail decay parameter at terminal time t (pt_plus > 0)
%       pt_minus - Left-tail decay parameter at terminal time t (pt_minus > 0)
%       ps_plus  - Right-tail decay parameter at initial time s (ps_plus > 0)
%       ps_minus - Left-tail decay parameter at initial time s (ps_minus > 0)
%
%   Outputs:
%       lambda   - Total expected number of jumps in the interval [s, t]
%
%   References:
%       Based on the finite-activity increment distribution framework (Baviera, 2026).

    % Reciprocal of the probability atom size (1/c).
    inverse_atom = (ps_plus .* ps_minus) ./ (pt_plus .* pt_minus);

    % lambda = -log(c) = log(1/c).
    lambda = log(inverse_atom);
end