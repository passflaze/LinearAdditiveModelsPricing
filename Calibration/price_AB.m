function CallPrice = price_AB(params, discount_factor, yf, sigma_ATM, moneyness_modified)
%PRICE_AB European call prices under the Additive Bachelier (AB) model.
%
% Vectorized, matrix-in / matrix-out pricing engine, built to consume the
% surface produced by moneyness_generator exactly like price_MA / price_GL.
% There is NO per-maturity for-loop: the normalized AB call price
% G(chi; eta, k) (Baviera & Massaria 2026, paper 3 Eq. 19) is
% maturity-independent under separability, so the whole (M x N) surface is
% obtained from a single Lewis-FFT evaluated on the full moneyness matrix,
% and the dollar price is recovered row-wise via Eq. (20):
%
%       C^mod(K,t) = B(t) * sigma_ATM(t) * sqrt(t) * G(chi; eta, k)
%
% INPUTS:
%   params             : 2-element vector [k, eta]
%                          k   > 0  : Inverse-Gaussian subordinator parameter
%                          eta      : asymmetry/skew parameter
%   discount_factor    : (M x 1) discount factors per maturity
%   yf                 : (M x 1) year fractions per maturity
%   sigma_ATM          : (M x 1) Bachelier ATM implied vols per maturity
%   moneyness_modified : (M x N) normalized moneyness chi = (K-F)/(sigma_ATM*sqrt(t)),
%                        with NaN where no market quote is selected.
%
% OUTPUT:
%   CallPrice          : (M x N) model call prices on the same grid as the
%                        input moneyness (NaN where moneyness is NaN).

    % 1. Extract parameters
    k   = params(1);
    eta = params(2);
    I_0=I0(0,k,eta); % I0 is the scale factor for the AB model, used to recover sigma_t from sigma_ATM via sigma_t = sigma_ATM / I_0.

    % 2. Normalized model price G(chi; eta, k) for the WHOLE surface.
    %    call_AB_FFT preserves the shape of its input and returns NaN where
    %    the moneyness is NaN, so a single FFT prices the (M x N) matrix.
    G = call_AB_FFT(moneyness_modified, k, eta, I_0);

    % 3. Maturity-dependent factor B * sigma_ATM * sqrt(t)  (M x 1), Eq. (20).
    maturity_multiplier = discount_factor .* sigma_ATM .* sqrt(yf);

    % 4. Dollar call surface via implicit expansion (M x 1) .* (M x N).
    CallPrice = maturity_multiplier .* G;

end
