function verbose_print(flag, varargin)
% Print only when flag is true.
if flag && ~isempty(varargin)
    fprintf(varargin{:});
end
end