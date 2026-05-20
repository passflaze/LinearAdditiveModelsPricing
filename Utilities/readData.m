function [strikes, calls, puts, expiries] = readData(callPath, putPath, snapDate, expiryFile)
% READDATA  Load WTI call and put option prices for one snapshot date.
%
%   [strikes, calls, puts, expiries] = READDATA(callPath, putPath, snapDate)
%   reads every option CSV in callPath and putPath and extracts the price
%   row corresponding to the value date snapDate. Each file holds the
%   option chain of one maturity; rows are returned sorted by expiry.
%
%   [...] = READDATA(..., expiryFile) also matches every file to its
%   expiry date, read from Expiries_Futures.txt via READEXPIRIES. The CSV
%   file name is the contract year-month (e.g. 202006.csv).
%   When the expiry is known, maturities that have already expired at
%   snapDate (expiry <= snapDate) are discarded, since they carry no
%   information for the calibration at the value date.
%
%   Untraded options keep their NaN price: the penny / no-arbitrage
%   filtering is left to the bootstrap stage.
%
%   Inputs
%     callPath   : folder containing the call CSV files
%     putPath    : folder containing the put CSV files
%     snapDate   : value date to extract, as datetime or yyyymmdd double
%     expiryFile : (optional) path to Expiries_Futures.txt
%
%   Outputs
%     strikes    : 1-by-nK strike grid (index points)
%     calls      : nT-by-nK call mid-prices, one row per maturity
%     puts       : nT-by-nK put mid-prices, aligned with calls
%     expiries   : nT-by-1 datetime expiry of each row (NaT if expiryFile
%                  is not provided)

% normalise the snapshot date: integer key for the CSV row lookup
% and a datetime for the expiry comparison
if isdatetime(snapDate)
    snapDay = dateshift(snapDate, 'start', 'day');
else
    snapDay = datetime(num2str(snapDate), 'InputFormat', 'yyyyMMdd');
end
snapKey = year(snapDay)*1e4 + month(snapDay)*1e2 + day(snapDay);

%list the option files 
callFiles = dir(fullfile(callPath, "*.csv"));
if isempty(callFiles)
    error("readData:noFiles", "no CSV files found in %s", callPath);
end
nT = numel(callFiles);

% match each file to its expiry date (just for clarity and to drop expired maturities, if any) 
expiries = NaT(nT, 1);
if nargin >= 4
    expTable = readExpiries(expiryFile);
    for k = 1:nT
        [~, name] = fileparts(callFiles(k).name);
        idx = find(expTable.yyyymm == str2double(name), 1);
        if isempty(idx)
            error("readData:noExpiry", ...
                  "no expiry matches file %s", callFiles(k).name);
        end
        expiries(k) = expTable.expiry(idx);
    end
end

% read every file and pick the snapshot row 
strikes = [];
calls   = [];
puts    = [];
for k = 1:nT
    cFile = fullfile(callFiles(k).folder, callFiles(k).name);
    pFile = fullfile(putPath, callFiles(k).name);   % same name in putPath
    if ~isfile(pFile)
        error("readData:missingPut", ...
              "no put file matching %s", callFiles(k).name);
    end

    Mc = readmatrix(cFile);
    Mp = readmatrix(pFile);

    % first row is the strike grid (with a leading NaN in column 1)
    if k == 1
        strikes = Mc(1, 2:end);
        calls   = nan(nT, numel(strikes));
        puts    = nan(nT, numel(strikes));
    elseif ~isequal(Mc(1, 2:end), strikes)
        error("readData:strikeMismatch", ...
              "%s has a different strike grid", callFiles(k).name);
    end

    % data rows: column 1 holds the value date as yyyymmdd
    row = find(Mc(2:end, 1) == snapKey, 1);
    if isempty(row)
        error("readData:snapNotFound", ...
              "value date %d not found in %s", snapKey, callFiles(k).name);
    end

    calls(k, :) = Mc(row + 1, 2:end);
    puts(k, :)  = Mp(row + 1, 2:end);
end

% sort by expiry and drop expired maturities (if expiries known)
if nargin >= 4
    [expiries, order] = sort(expiries);
    calls = calls(order, :);
    puts  = puts(order, :);

    alive    = expiries > snapDay;     % keep only not-yet-expired maturities
    expiries = expiries(alive);
    calls    = calls(alive, :);
    puts     = puts(alive, :);
end
end
