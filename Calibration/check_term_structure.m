function [ok, bad_idx] = check_term_structure(sigma_ATM, yf, expiries)
% CHECK_TERM_STRUCTURE  Verify sigma_ATM*sqrt(t) is strictly increasing in t,
%   the additive-model precondition (a data condition no (alpha,beta) can fix).
%
% INPUTS
%   sigma_ATM, yf : (M x 1) ATM vols and year fractions (sorted by maturity)
%   expiries      : (optional) (M x 1) datetime, for warning labels
% OUTPUTS
%   ok      : true iff strictly increasing
%   bad_idx : indices (>=2) where monotonicity fails (empty if ok)

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
