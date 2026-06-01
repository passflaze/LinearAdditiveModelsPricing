function price = pricing_fwd_start_GL_analytic(params, sc_T1, sc_T2, df, K2, ...
                                               forward, fwd_factor, M, dz)
%PRICING_FWD_START_GL_ANALYTIC  Semi-analytic forward-start price under the
%   Generalized Logistic (GL) model.
%
%   Prices  df * E[ ( S(T2) - K2*F(T1,T2) )_+ ]  by Lewis-FFT inversion of the
%   Lemma-2 (Forward.pdf) increment CF cf_increment_GL, conditioning on the
%   reset forward. See fwd_start_LA_analytic for the methodology. At K2 = 1 it
%   returns the ATM forward-start df * E[max(W,0)] (matches the Step-8 uniform
%   Lewis call in run_ex3).
%
% INPUTS
%   params     : GL parameters [alpha; beta]
%   sc_T1,sc_T2: FULL scale factors (sigma_ATM/I_0)*sqrt(T) at T1 and T2
%   df         : B(0,T2) discount factor
%   K2         : strike multiplier(s), scalar or vector
%   forward    : F(0,T2)
%   fwd_factor : Lemma 2 rescaling B(0,T1)/B(0,T2)  (optional, default 1)
%   M, dz      : Lewis-FFT grid parameters             (optional, default 16, 0.05)
%
% OUTPUT
%   price      : 1 x numel(K2) row of forward-start prices

    if nargin < 7, fwd_factor = []; end
    if nargin < 8, M  = []; end
    if nargin < 9, dz = []; end

    price = fwd_start_LA_analytic(@cf_increment_GL, 'GL', params, sc_T1, sc_T2, ...
                                  df, K2, forward, fwd_factor, M, dz);
end
