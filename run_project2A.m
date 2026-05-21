% data loading
clear,clc
addpath("Utilities/");

callpath="Data/datacalls";
putpath="Data/dataputs";
expiryFile = "Data/Expiries_Futures.txt";
valueDate = datetime(2020,06,02);

[strikes, calls, puts, expiries] = readData(callpath, putpath, valueDate, expiryFile);

%% bootstrap: synthetic discount factor and forward, one per maturity
nT = numel(expiries);
discount_factor = zeros(nT,1);
forward = zeros(nT,1);
R2 = zeros(nT,1);
call_atm = zeros(nT,1);

for k = 1:nT
    [discount_factor(k), forward(k), R2(k)] = bootstrap(puts(k,:), calls(k,:), strikes);
    call_atm(k) = callATM(calls(k,:), puts(k,:), strikes, forward(k), discount_factor(k));
end

act_365 = 3; % actual/365 day count convention for time to maturity
yf = yearfrac(valueDate, expiries, act_365); 
sigma_atm = sigmaATM(call_atm, discount_factor, yf,expiries);

