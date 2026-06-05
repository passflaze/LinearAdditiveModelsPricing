function [c, ceq] = term_structure_nonlcon(~, sigma_ATM, yf)
% TERM_STRUCTURE_NONLCON  fmincon constraint: sigma_t*sqrt(t) increasing in t
%   (additive-model precondition; data-only, independent of the parameters).
%
% INPUTS:
%   ~          : parameter vector (unused; required by the fmincon nonlcon signature)
%   sigma_ATM  : (M x 1) ATM vols per maturity
%   yf         : (M x 1) year fractions per maturity
% OUTPUTS:
%   c          : (M-1 x 1) inequality constraints, c(i) = ts(i) - ts(i+1) <= 0
%                enforce ts = sigma_ATM*sqrt(t) strictly increasing
%   ceq        : empty (no equality constraints)

    ts  = sigma_ATM(:) .* sqrt(yf(:));
    c   = ts(1:end-1) - ts(2:end);   % <= 0 enforces increasing
    ceq = [];

end
