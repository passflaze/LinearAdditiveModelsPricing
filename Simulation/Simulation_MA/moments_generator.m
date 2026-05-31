function [mean_val, var_val, skew_val, kurt_val] = moments_generator(pt_plus, pt_minus, ps_plus, ps_minus, drift)
% MOMENTS_GENERATOR Computes the exact analytical first 4 moments of the 
% Minimal Additive (MA) increment distribution using cumulants.
%
% Inputs:
%   pt_plus, pt_minus : Right and left tail decay parameters at time t
%   ps_plus, ps_minus : Right and left tail decay parameters at time s
%   drift             : Deterministic drift component (Delta_mu)
%
% Outputs:
%   mean_val : Theoretical Mean (First Moment)
%   var_val  : Theoretical Variance (Second Central Moment)
%   skew_val : Theoretical Skewness (Third Standardized Moment)
%   kurt_val : Theoretical Pearson Kurtosis (Fourth Standardized Moment)

    % Cumulant 1 (Mean)
    k1 = drift + (1./pt_plus) - (1./pt_minus) - (1./ps_plus) + (1./ps_minus);
    
    % Cumulant 2 (Variance)
    k2 = (1./(pt_plus.^2)) + (1./(pt_minus.^2)) - (1./(ps_plus.^2)) - (1./(ps_minus.^2));
    
    % Cumulant 3
    k3 = 2 * ( (1./(pt_plus.^3)) - (1./(pt_minus.^3)) - (1./(ps_plus.^3)) + (1./(ps_minus.^3)) );
    
    % Cumulant 4
    k4 = 6 * ( (1./(pt_plus.^4)) + (1./(pt_minus.^4)) - (1./(ps_plus.^4)) - (1./(ps_minus.^4)) );

    % Map Cumulants to Standard Statistical Moments
    mean_val = k1;
    var_val  = k2;
    skew_val = k3 ./ (k2 .^ 1.5);
    kurt_val = (k4 ./ (k2 .^ 2)) + 3; % Standard kurtosis (Normal = 3)
end