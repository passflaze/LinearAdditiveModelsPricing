function varargout = nth_out(idx, fcn, varargin)
%NTH_OUT  Return a selected output of a multi-output function.
%   y = NTH_OUT(N, FCN, ARGS...) calls FCN(ARGS...) requesting N outputs and
%   returns only the N-th. Useful to grab a deep output (e.g. the MC sample
%   std at output position 5) of a function inside an anonymous-function
%   closure, where the usual [~,~,~,~,s] = fcn(...) syntax is not available.
%
%   Example:
%       sig = @(N) nth_out(5, @CoC_pricing_MC, params, scale, N, ...);
%       s   = sig(1e4);   % runs the pilot and returns only sigma

    out = cell(1, idx);
    [out{:}] = fcn(varargin{:});
    varargout = out(idx);
end
