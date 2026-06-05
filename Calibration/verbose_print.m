function verbose_print(flag, varargin)
% VERBOSE_PRINT  Conditional fprintf wrapper.
%
% INPUTS:
%   flag     : logical; print only when true
%   varargin : fprintf-style format string and arguments
% OUTPUT:
%   none (writes to the command window)
if flag && ~isempty(varargin)
    fprintf(varargin{:});
end
end