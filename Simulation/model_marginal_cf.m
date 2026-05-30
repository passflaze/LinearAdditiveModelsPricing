function [phi, c_neg, c_pos] = model_marginal_cf(model, params, T, sigma_T)
% MODEL_MARGINAL_CF  Marginal characteristic function and analyticity strip
% of a Linear Additive log-price (forward increment) at maturity T.
%
% Single dispatch point for the AB / GL simulation pipeline: every other
% function in Simulation/ is model-agnostic and only depends on the handle
% and the two strip constants returned here.
%
% INPUTS
%   model    : 'AB' (Additive Bachelier) or 'GL' (Generalized Logistic)
%   params   : struct with the model parameters
%                AB -> params.k , params.eta
%                GL -> params.alpha , params.beta
%   T        : maturity / reset date (year fraction). T = 0 -> degenerate CF = 1.
%   sigma_T  : calibrated scale at T (= sigma_ATM(T)/I_0).
%
% OUTPUTS
%   phi      : handle @(u) -> marginal CF of f_{T,T}, vectorised over u
%              (u may be complex, for the Lewis contour shift).
%   c_neg    : strip half-decay controlling the NEGATIVE contour shift a_neg
%              (right-tail accuracy). Admissible shift = c_neg/(sigma_T*sqrt(T)).
%   c_pos    : strip half-decay controlling the POSITIVE contour shift a_pos
%              (left-tail accuracy).  Admissible shift = c_pos/(sigma_T*sqrt(T)).
%
% References: AB strip p_+- = -+eta + sqrt(eta^2 + 1/k) (ex2/I0_AB, Baviera-
% Manzoni 2026 Sec.5.3). GL strip u in (-alpha/sT, beta/sT) (ex3_GL/lewis_fft_cdf).

    switch upper(model)
        case 'AB'
            k   = params.k;
            eta = params.eta;
            phi = charateristic_function_AB(T, k, eta, sigma_T);   % returns @(u) col
            c_neg = eta + sqrt(eta^2 + 1/k);    % p_+  -> a_neg
            c_pos = -eta + sqrt(eta^2 + 1/k);   % p_-  -> a_pos

        case 'GL'
            alpha = params.alpha;
            beta  = params.beta;
            sT    = sigma_T * sqrt(T);
            % cf_GL signature (Distributions/cf_GL.m): cf_GL(u, [alpha,beta], scale)
            % -> pass the param VECTOR and let cf_GL apply the scale internally
            % (u is scaled by sT inside). Same convention as cf_increment_GL.
            phi   = @(u) cf_GL(u, [alpha, beta], sT);
            c_neg = alpha;   % left edge of (-alpha, beta) -> a_neg
            c_pos = beta;    % right edge                  -> a_pos

        otherwise
            error('model_marginal_cf:unknownModel', ...
                  'Unknown model "%s". Use ''AB'' or ''GL''.', model);
    end
end
