function [ok, bad_idx] = check_term_structure(sigma_ATM, yf, expiries)
% CHECK_TERM_STRUCTURE Verify the additive-model precondition on sigma_ATM.
%
% The MA / GL family (cf. MA_Basic_Project, Sec. 1.3 and paper 3 Sec. 2)
% requires sigma_t * sqrt(t) to be (strictly) increasing in t so that the
% rescaled tails p_t^+- = (alpha or beta)/(sigma_t*sqrt(t)) are decreasing.
% With (alpha, beta) constant across maturities, sigma_t = sigma_ATM_t / I0
% and I0 does not depend on t, so the requirement reduces to
%
%       sigma_ATM_t * sqrt(t)   strictly increasing in t.
%
% This is a *data* condition: no choice of (alpha, beta) can fix it.
%
% INPUTS
%   sigma_ATM : (M x 1) Bachelier ATM vols, one per maturity
%   yf        : (M x 1) year fractions, sorted by maturity
%   expiries  : (optional) (M x 1) datetime vector, used only for warning
%               labels
%
% OUTPUTS
%   ok      : true iff sigma_ATM .* sqrt(yf) is strictly increasing
%   bad_idx : indices i (>=2) where the monotonicity fails
%             (i.e. ts(i) <= ts(i-1)); empty when ok is true

    sigma_ATM = sigma_ATM(:);
    yf        = yf(:);

    ts  = sigma_ATM .* sqrt(yf);
    dts = diff(ts);
    bad_idx = find(dts <= 0) + 1;
    ok = isempty(bad_idx);

    if ok
        fprintf(['  -> Term structure check: sigma_ATM*sqrt(t) strictly ' ...
                 'increasing across %d maturities. OK.\n'], numel(ts));
        return
    end

    if nargin >= 3 && ~isempty(expiries)
        labels = strjoin(cellstr(string(expiries(bad_idx), 'yyyy-MM-dd')), ', ');
    else
        labels = strjoin(arrayfun(@(k) sprintf('t=%g', yf(k)), bad_idx, ...
                                  'UniformOutput', false), ', ');
    end

    warning('AddBach:TermStructure', ...
        ['sigma_ATM*sqrt(t) is not strictly increasing at maturities: %s. ' ...
         'The additive-model constraint p_t^+- decreasing is violated. ' ...
         'No (alpha, beta) calibration can repair this -- it is a data condition.'], ...
        labels);

end
