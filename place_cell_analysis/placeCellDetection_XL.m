clear all 
close all 
addpath C:\Users\xl313\OneDrive\Documents\GitHub\spike-analysis\place_cell_analysis
addpath('C:\Users\xl313\OneDrive\Documents\GitHub\spike-analysis\npy-matlab-master\npy-matlab')
addpath('C:\Users\xl313\OneDrive\Documents\GitHub\spike-analysis\extractWF_xl')

%% Load smoothed keypoint predictions
behavior_rootDir = 'Z:\Sherry\poseTrackingXL\training_files\raw_acquisition_copy\RBY52_012025';
% to load in x,y positions of COM body centorid (xCom, yCom)
load(fullfile(behavior_rootDir,'behavioral_data\xyzCom_022025.mat')); % "com_pos_smooth","com_vel_smooth"
% verify that frames across all cameras are in sync
addpath('Z:\Sherry\camera');
sync = verifyBehavioralVideoMetadata(behavior_rootDir); % sync contains info from each video csv file

%% Extract timestamps for xy positions for the selected period
tt = sync{2,1}(1:length(com_pos_smooth),2);
tt = (tt-tt(1))/10e8; % now in sec, IFI = 0.02
IFI = diff(tt);
fprintf("Overall video acqusition frame rate: %.3f\n", 1/mean(IFI));
% variability of IFI
ecdf(IFI);

%%
B = struct(...
    'body_xpos', com_pos_smooth(:,2,1), ...
    'body_ypos', com_pos_smooth(:,2,2), ...
    'body_zpos', com_pos_smooth(:,2,3), ...
    'tt', tt' ...s
);
fps = double(1/mean(diff(B.tt)));
% Compute velocity (mm/sec) manually [(position_t - position_t-1) / Delta_t]
B.speed = [sqrt(diff(B.body_xpos).^2+diff(B.body_ypos).^2 + diff(B.body_zpos).^2)*fps; NaN];

%% Visualisation: behavioral information 
% To plot the x, y, z coordinates of the body centroid of the bird over frames
figure;
hold on 
plot(B.body_xpos, 'LineWidth', 1.2, 'DisplayName', 'x', 'LineStyle', '-');
plot(B.body_ypos, 'LineWidth', 1.2, 'DisplayName', 'y', 'LineStyle', '-');
plot(B.body_zpos, 'LineWidth', 1.2, 'DisplayName', 'z', 'LineStyle', '-');
plot(B.speed, 'LineWidth', 1.2, 'DisplayName', 'V (mm/sec)', 'LineStyle', '-');
hold off
legend('Location', 'northeastoutside');
title('Body Centroid (ComNet)');
xlabel('Frames (50Hz)')
ylabel('Position (mm)')
yline(300, '--', 'Color', 'r', 'Label', 'y = 300', DisplayName = 'Spd Limit');

%% Mask low speed (based on Hannah's processSpeed.m script)
speed_thresh = 300; %mm/s
max_duration = 5; % (sec)
code_area = regionprops(B.speed<speed_thresh,'area'); % duration (frames) of < speed thresh 
code_label = bwlabel(B.speed<speed_thresh);
stationary_lengths = [code_area.Area]/fps; % the corresponding duration 
reject_segs = find(stationary_lengths > max_duration); % stationary for more than certain time
B.speed_mask = false(size(B.tt));
for jj = 1:length(reject_segs)
    B.speed_mask(code_label==reject_segs(jj)) = 1;
end

figure; plot(B.tt,B.speed,'k'); hold on
temp = B.speed; temp(B.speed_mask) = NaN;
plot(B.tt, temp,'r');   xlabel('Time (s)'); ylabel('Speed (mm/s)')
fprintf('%f excluded by speed thresholding\n', mean(B.speed_mask))

% Apply the speed thresholding to the output
B.speed_raw = B.speed; % Keep a copy to analyze speed tuning
B.body_xpos(B.speed_mask) = NaN;
B.body_ypos(B.speed_mask) = NaN;
B.body_zpos(B.speed_mask) = NaN;
B.speed(B.speed_mask) = NaN;

%% Visualise smoothed data with and without speed thresholding 
% Extract the unmasked x and y coordinates
xPos0 = com_pos_smooth(:, 2, 1); 
yPos0 = com_pos_smooth(:, 2, 2); 

num_frames = size(com_pos_smooth, 1); % Number of frames
colors = linspace(0, num_frames, num_frames); 
colormap('parula'); 
set(gcf,'Position', [100, 100, 1200, 500]); % [left, bottom, width, height]

% Plot the 2D positions
subplot(1,2,1);
scatter(xPos0, yPos0, 20, colors, 'filled', 'MarkerFaceAlpha', 0.7); % Scatter plot
xlabel('X (mm)');
ylabel('Y (mm)');
xlim([-350, 350]); 
ylim([-350, 350]); 
title('Bird in Arena (Smoothed)');
grid on;

% Visualise occupancy with speed > 5cm/s
num_frames = size(com_pos_smooth, 1);
colors = linspace(0, num_frames, num_frames);
colormap('parula'); 

% Plot the 2D positions
subplot(1,2,2);
scatter(B.body_xpos, B.body_ypos, 20, colors, 'filled', 'MarkerFaceAlpha', 0.7); % Scatter plot
xlabel('X (mm)');
ylabel('Y (mm)');
xlim([-350, 350]); 
ylim([-350, 350]); 
title('Bird in Arena (V > 300mm/s)');
grid on;

colorbar;
caxis([0, num_frames]); 
colorbar('Ticks', [0, num_frames/2, num_frames], 'TickLabels', {'Start', 'Mid', 'End'}); % Customize colorbar ticks
ylabel(colorbar, 'Frames'); 

%% Read in Ephys file  
filepath = 'Z:\Sherry\ephys_acquisition\RBY52\4turnsTotal2ndPart_250120_122553\';

% 0. Check synchronisation between cameras and Intan recording (TTL pulses)
% compute the TTL pulses recorded in Intan file (dig_time for the entire duration)
[dig_raw_all, h_all, dig_time_all] = readIntanDig_In1only(filepath);
TTL_all = cellfun(@length, dig_time_all);
fprintf("The total number of TTL pulses for the entire recording: %d\n", TTL_all)
frames_all =cellfun(@length, sync); 
fprintf("The number of frames in 1 bottom and 4 side cameras: %d\n", frames_all)
if ~all(frames_all(2:end)==TTL_all)
    warning("different TTL counts between cameras and ephys")
end
% check if there is any IFI of TTL pulses larger than 1ms (0.001s), if so,
% find the indices
target_TTL_rate = 1/50;
if any(abs(diff(dig_time_all{1,1}(:,1)) - target_TTL_rate) > 0.001) % dig_time_all{1,1} is rising phase
    long_idx = find(abs(diff(dig_time_all{1,1}(:,1)) - target_TTL_rate) > 0.001); 
    dig_time_IFI = diff(dig_time_all{1,1}(:,1));
    long_IFI = dig_time_IFI(long_idx);
    fprintf("Indices and values of inter TTL pulses intervals > 1ms: %d, %.4fs\n",long_idx, long_IFI);
end

%% 1. Read in 1hr of Intan digitalIn.dat file
tstart = 0;
tduration = 3600;  % in seconds, % so to extract 1 hr, t = 3600 
[dig_raw, h, dig_time] = readIntanDig_In1only(filepath, tduration, tstart);
% dig_raw : Matrix containing the raw digital data.(n channels x n samples).
% h       : Struct containing header information from the info.rhd file.
% dig_time: Cell array, each element corresponds to one digitalinput, and contains start (1st col) and stop (2nd col) 
%           timestamps for the digital input
TTL_rate = length(dig_time{1})/(tduration-dig_time{1,1}(1));
fprintf("Rate of TTL pulses: %.3f\n", TTL_rate); % no. of TTL pulses/TTL_duration

%% 2. Read in Kilosort data 
ksDir = 'Z:\Sherry\ephys_acquisition\RBY52\4turnsTotal2ndPart_250120_122553\kilosort4';
fs = 30000; % sample rate
option_only_good = 1;
session_label = '4turnsTotal2ndPart_250120_122553';

% Read in spike times 
st = readNPY(fullfile(ksDir,'spike_times.npy')); % in samples 
st = double(st)/fs; % in seconds 
% Read in unit IDs for each spike, and extract "good" units 
sID = readNPY(fullfile(ksDir,'spike_clusters.npy'));
[unit_ID, goodLabels] = getPhyClusterLabels(ksDir); 
fprintf("Numbers of good & MUA units: %d\n Number of good units only: %d\n", length(goodLabels), sum(strcmp(goodLabels, 'good')));

% Finding the indices of the "good" units only 
inds_keep = find(strcmp(goodLabels,'good'));

%%  3. Loop through units.
S = struct;
ii_out = 1;
for ii = inds_keep(:)' 
    % Cell ID
    cID = unit_ID(ii); % select out the "good" unit IDs
    % Get spike times from this unit
    tt = st(sID==cID); 
    % Store spikes
    tstart = min(tt);
    tend = max(tt);

    % Store this cell
    S(ii_out).data = tt;
    S(ii_out).tstart = tstart;
    S(ii_out).tend = tend;
    S(ii_out).units = 's';
    S(ii_out).cellID = cID;
    S(ii_out).TTL_rate = TTL_rate;
    ii_out = ii_out+1;
end 

%% 4. extract spikes in a defined period 
TTL_start = dig_time{1,1}(1,1); % the timestamp when TTL pulses start
TTL_end = dig_time{1,1}(end,1) + TTL_start; 
% to make TTL timestamps start at 0 to match frame timestamps
TTL_0start = dig_time{1,1}(:,1)- TTL_start;

% find the corresponding frame indices for each spike 
st_unit_onlyTTL = cell(1, size(S,2));
ISPK_idx = cell(1,size(S,2));  
for i = 1:size(S,2) % the number of good units
    % to only include the spikes happening within TTL start-end period
    st_unit_onlyTTL{i} = S(i).data(S(i).data > TTL_start & S(i).data < TTL_end);
    % since counting starts at the first TTL pulse, so subtract the first TTL timestamp from spike times 
    st_unit_onlyTTL{1,i} = st_unit_onlyTTL{1,i}-TTL_start;
    % find the indices of TTLs (frames) where the spikes happen 
    indices = NaN(1, length(st_unit_onlyTTL{1,i}));
    for j = 1:length(st_unit_onlyTTL{1,i})
          indices(1,j)= find(TTL_0start <= st_unit_onlyTTL{1,i}(j),1,'last');
    end
    ISPK_idx{i} = indices;
end

%% Visualise each unit's spikes in a raster plot
figure;
hold on;
TTL_rise = dig_time{1,1}(:,1);

for trial = 1:length(st_unit_onlyTTL)
    scatter(st_unit_onlyTTL{trial}, trial * ones(size(st_unit_onlyTTL{trial})), 8, 'k', 'filled');
end
xlabel('Time (s)');
ylabel('Neurons')
hold off;

%% 5. placeCellAnalysis.m  
%   XX and YY: animal's position (mm)
xx = B.body_xpos;
yy = B.body_ypos;
FPS = S.TTL_rate; %   FPS: frames per second of the behavior % In reality it is not actually 50!

savedir = fullfile(behavior_rootDir, 'place_maps_v2_noLowSpeed_aligned');
if ~exist(savedir, 'dir')
    mkdir(savedir);
end 

for i = 1:size(st_unit_onlyTTL,2) 
    ISPK = ISPK_idx{1,i};  %   ISPK: spike indices into the behavior ISPK (10, 23, 37)
    NSHUFFLE = 0; %   NSHUFFLE: number of random shuffles to conduct (specify 0 for no shuffling)
    LIMS = [-330.2 330.2]; %   LIMS: two-number vector limits of the square arena (*[0 609.6] (mm) in the 2-foot arena)
    NBINS = 40; %   NBINS: number of bins in which to divide the behavior *(40 bins)
    NHAMMING = 13; %   NHAMMING: size of hamming filter *(13)
    MIN_OCCUP = 0.1; %   MIN_OCCUP: minimum occupancy time in each bin *(0.1 s)
    [map, info, bin_centers, pos_time_smooth] = placeCellAnalysis(xx, yy, FPS, ISPK, NSHUFFLE, LIMS, NBINS, NHAMMING, MIN_OCCUP,1);
    title(sprintf('Cell ID: %d',S(i).cellID));
    colorbar;
    filename = fullfile(savedir,sprintf('cell%d.png', i));
    saveas(gcf,filename);
end 

%%
lims = [-330.2 330.2];
bin_edges = linspace(lims(1), lims(2), 40+1);
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;

% Calculate position counts (doesn't change with shuffling)
pos_count = histcounts2(xx, yy, bin_edges, bin_edges)'; %% 
xx_sp = xx(ISPK);
yy_sp = yy(ISPK);
sp_count = histcounts2(xx_sp,yy_sp, bin_edges, bin_edges)';

% Filter while dealing with NaNs (file exchange function nanconv)
sp_count_smooth = nanconv(sp_count, w,'nanout', 'edge');
curr_map = sp_count_smooth./pos_time_smooth;

