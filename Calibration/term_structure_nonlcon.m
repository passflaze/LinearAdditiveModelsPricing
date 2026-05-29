function [c, ceq] = term_structure_nonlcon(~, sigma_ATM, yf)
% TERM_STRUCTURE_NONLCON  fmincon nonlinear constraint encoding the
% additive-model preconditions (cf. MA_Basic_Project, Sec. 1.3):
%
%       sigma_t * sqrt(t)  strictly increasing in t
%       p_t^+- = (alpha,beta)/(sigma_t*sqrt(t))  decreasing in t
%
% With a single (alpha, beta) pair calibrated across all maturities the
% second condition follows from the first (alpha, beta > 0 constant
% implies p^+- decreasing iff sigma_t*sqrt(t) increasing). The first
% condition only depends on sigma_ATM, not on the parameters; we still
% encode it as `nonlcon` so that:
%
%   (a) the call to fmincon documents the modeling assumption, and
%   (b) extending the calibration to slice-by-slice (alpha_t, beta_t)
%       lights up a constraint that actually binds, without changing the
%       run script.
%
% fmincon nonlcon convention:
%   c(x)   <= 0   (inequality)
%   ceq(x) == 0   (equality)
%
% INPUTS:
%   x         : (unused, first arg required by fmincon nonlcon signature)
%   sigma_ATM : (M x 1) Bachelier ATM vols
%   yf        : (M x 1) year fractions
%
% OUTPUTS:
%   c   : (M-1) x 1 vector. c(i) = ts(i) - ts(i+1)  must be <= 0,
%         i.e. ts strictly increasing.
%   ceq : empty.

    ts  = sigma_ATM(:) .* sqrt(yf(:));
    c   = ts(1:end-1) - ts(2:end);   % <= 0 enforces increasing
    ceq = [];

end
