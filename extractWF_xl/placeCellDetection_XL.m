addpath D:\GitHub\spike-analysis\extractWF_xl

%% Read in Ephys file  
filepath = 'Z:\Sherry\ephys_acquisition\RBY52\4turnsTotal2ndPart_250120_122553\';
tduration = 1800;  % in seconds, % so to extract 0.5 hr, t = 1800 
tstart = 0;

%% Read in Intan digitalIn.dat file
[dig_raw, h, dig_time] = readIntanDig_In1only(filepath, tduration,tstart);

fprintf("Acquisition frame rate (FRS): %.3f\n", length(dig_time{1})/tduration);


%% Specify parameters for placeCellAnalysis function 
xx = com_pos_smooth(:,1,1); %(:,:,1) is the body centroid
yy = com_pos_smooth(:,2,1);
%   XX and YY: animal's position (mm)
FPS = 50; %   FPS: frames per second of the behavior


ISPK = dig_time{1,1}; %   ISPK: spike indices into the behavior ISPK

NSHUFFLE = 0; %   NSHUFFLE: number of random shuffles to conduct (specify 0 for no shuffling)
LIMS = [0,660.4]; %   LIMS: two-number vector limits of the square arena (*[0 609.6] (mm) in the 2-foot arena)
NBINS = 40; %   NBINS: number of bins in which to divide the behavior *(40 bins)
NHAMMING = 13; %   NHAMMING: size of hamming filter *(13)
MIN_OCCUP = 0.1; %   MIN_OCCUP: minimum occupancy time in each bin *(0.1 s)


placeCellAnalysis(xx, yy, FPS, ISPK, NSHUFFLE, LIMS, NBINS, NHAMMING, MIN_OCCUP)