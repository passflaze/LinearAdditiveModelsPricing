function [c, ceq] = term_structure_nonlcon(~, sigma_ATM, yf)
% TERM_STRUCTURE_NONLCON  fmincon constraint: sigma_t*sqrt(t) increasing in t
%   (additive-model precondition; data-only, independent of the parameters).
%   c(i) = ts(i) - ts(i+1) <= 0 enforces the monotonicity; ceq empty.
%   First arg unused (required by the fmincon nonlcon signature).

    ts  = sigma_ATM(:) .* sqrt(yf(:));
    c   = ts(1:end-1) - ts(2:end);   % <= 0 enforces increasing
    ceq = [];

end
