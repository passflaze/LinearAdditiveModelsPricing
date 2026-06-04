function CallPrice = price_AB(params, discount_factor, yf, sigma_ATM, moneyness_modified)
%PRICE_AB  European call surface under the Additive Bachelier model.
%   By separability the normalized call G(chi; eta, k) is maturity-independent
%   ([3] Eq. 19), so one Lewis-FFT prices the whole (M x N) surface and the
%   dollar price follows row-wise ([3] Eq. 20):
%       C(K,t) = B(t) * sigma_ATM(t) * sqrt(t) * G(chi).
%
% INPUTS:
%   params             : [k, eta]  (k > 0 IG parameter, eta skew)
%   discount_factor, yf, sigma_ATM : (M x 1) per maturity
%   moneyness_modified : (M x N) chi = (K-F)/(sigma_ATM*sqrt(t)), NaN if no quote
% OUTPUT:
%   CallPrice          : (M x N) call prices (NaN where chi is NaN)

    k   = params(1);
    eta = params(2);
    I_0 = I0_AB(0, params);                       % sigma_t = sigma_ATM / I_0

    G = call_AB_FFT(moneyness_modified, k, eta, I_0);   % normalized call (M x N)
    maturity_multiplier = discount_factor .* sigma_ATM .* sqrt(yf);   % (M x 1)
    CallPrice = maturity_multiplier .* G;

end
