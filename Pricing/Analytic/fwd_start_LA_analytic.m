function price = fwd_start_LA_analytic(cf_inc, model, params, sc_T1, sc_T2, ...
                                       df, K2, forward, fwd_factor, M, dz)
%FWD_START_LA_ANALYTIC  Semi-analytic forward-start price for the Linear
%   Additive AB/GL models, via Lewis-FFT pricing of the Lemma-2 increment.
%
%   Shared engine behind pricing_fwd_start_AB_analytic / _GL_analytic. With the
%   Lemma 2 (Forward.pdf) reconstruction
%       F(T1,T2) = forward + fwd_factor*Z1 ,   Z1 = f_{T1,T1}
%       S(T2)    = forward + fwd_factor*Z1 + W ,   W = f_{T2,T2} - f_{T1,T2}
%   and W independent of Z1, the forward-start payoff factorises as
%       [ S(T2) - K2*F(T1,T2) ]_+ = [ W - K_star*(forward + fwd_factor*Z1) ]_+
%   with K_star = K2 - 1. Conditioning on the reset forward A = forward+fwd*Z1,
%       price = df * E_{Z1}[ c_W( K_star*A ) ] ,
%   where c_W(kappa) = E[max(W - kappa, 0)] is the increment call priced by
%   lewis_FFT_call and Z1 ~ reset marginal f_{T1,T1}. At K2 = 1 (K_star = 0)
%   this collapses to the ATM forward-start df * c_W(0).
%
% INPUTS
%   cf_inc     : increment CF handle, signature cf(u, params, scale, fwd_factor)
%                (cf_increment_AB / cf_increment_GL)
%   model      : 'AB' | 'GL'  (selects the Lewis damping shift)
%   params     : model parameter column vector (AB [k; eta], GL [alpha; beta])
%   sc_T1,sc_T2: FULL scale factors (sigma_ATM/I_0)*sqrt(T) at T1 and T2
%   df         : B(0,T2) discount factor (payoff received at T2)
%   K2         : strike multiplier(s), scalar or vector
%   forward    : F(0,T2)
%   fwd_factor : Lemma 2 rescaling B(0,T1)/B(0,T2)  (optional, default 1)
%   M, dz      : Lewis-FFT grid parameters             (optional, default 16, 0.05)
%
% OUTPUT
%   price      : 1 x numel(K2) row of forward-start prices

    if nargin < 9  || isempty(fwd_factor), fwd_factor = 1;    end
    if nargin < 10 || isempty(M),          M  = 16;           end
    if nargin < 11 || isempty(dz),         dz = 0.05;         end

    % Wrap the Lemma-2 increment CF to the lewis_FFT_* interface cf(u,params,sc).
    cf_h  = @(u, p, sc) cf_inc(u, p, sc, fwd_factor);
    scvec = [sc_T1; sc_T2];
    K2    = K2(:).';
    price = zeros(1, numel(K2));

    % Reset-forward marginal (scale_factor(1) = 0 collapses cf_inc to phi_T1).
    % Only needed for non-ATM strikes; a coarse node set is plenty for the
    % smooth 1-D outer integral.
    if any(K2 ~= 1)
        [cdf_z1, x_z1] = lewis_FFT_digital(cf_h, M, dz, params, [0; sc_T1], ...
                                           true, model, true, false);
        idx    = round(linspace(1, numel(x_z1), min(4000, numel(x_z1))));
        x_z1   = x_z1(idx);  cdf_z1 = cdf_z1(idx);
        w_pmf  = diff(cdf_z1(:));                          % sums to 1 (cdf 0 -> 1)
        x_mid  = 0.5*(x_z1(1:end-1) + x_z1(2:end));
        A_node = forward + fwd_factor * x_mid(:);          % reset forward per node
    end

    for j = 1:numel(K2)
        Kstar = K2(j) - 1;
        if Kstar == 0
            % ATM forward-start: payoff = max(W, 0) -> single increment call at 0.
            price(j) = df * lewis_FFT_call(cf_h, M, dz, params, scvec, 0, true, model);
        else
            % General strike: integrate the increment call over the reset forward.
            kappa    = Kstar * A_node;                     % dollar strikes
            cW       = lewis_FFT_call(cf_h, M, dz, params, scvec, kappa, true, model);
            price(j) = df * sum(cW(:) .* w_pmf);
        end
    end
end
