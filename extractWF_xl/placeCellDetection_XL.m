addpath D:\GitHub\spike-analysis\extractWF_xl

% Specify parameters 
filepath = 'Z:\Sherry\ephys_acquisition\RBY52\4turnsTotal2ndPart_250120_122553\';
tduration = 3600;  % in seconds, % so to extract 1 hr, t = 3600 
tstart = 0;

% To read in Intan digitalIn.dat file
[dig_raw, h, dig_time] = readIntanDig_In1only(filepath, tduration,tstart);

fprintf("Acquisition frame rate (FRS): %.3f\n", length(dig_time{1})/tduration);

