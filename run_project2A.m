% data loading
addpath("Utilities/");

callpath="Data/datacalls";
putpath="Data/dataputs";
expiryFile = "Data/Expiries_Futures.txt";
today = datetime(2020,06,02);

[strikes, calls, puts, expiries] = readData(callpath, putpath, today, expiryFile);


