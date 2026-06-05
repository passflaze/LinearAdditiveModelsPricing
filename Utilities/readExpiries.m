function expTable = readExpiries(filePath)
% READEXPIRIES  Read the futures expiry file.
%
%   expTable = READEXPIRIES(filePath) parses Expiries_Futures.txt, whose
%   rows have the form "<CODE> <YYYY/MM/DD>" (e.g. "JUN20 2020/05/14").
%
% INPUTS:
%   filePath - path to Expiries_Futures.txt
%
% OUTPUT:
%   expTable - table with columns:
%                code   (string)   contract code, e.g. "JUN20"
%                expiry (datetime) expiry date
%                yyyymm (double)   contract year-month, e.g. 202006

raw = strip(readlines(filePath));
raw = raw(raw ~= "");                       % drop empty lines

parts  = split(raw, whitespacePattern);     % N-by-2 string array
code   = parts(:,1);
expiry = datetime(parts(:,2), 'InputFormat', 'yyyy/MM/dd');

% map the 3-letter month code to a month number
monthCodes = ["JAN" "FEB" "MAR" "APR" "MAY" "JUN" ...
              "JUL" "AUG" "SEP" "OCT" "NOV" "DEC"];
[found, mon] = ismember(extractBefore(code, 4), monthCodes);
if ~all(found)
    error("readExpiries:badCode", ...
          "unrecognised month code in: %s", strjoin(code(~found), ", "));
end

yr     = 2000 + double(extractAfter(code, 3));
yyyymm = yr*100 + mon;

expTable = table(code, expiry, yyyymm);
end
