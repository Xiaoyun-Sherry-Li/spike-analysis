function WFstruct = extractWaveforms(dataDir,ksDir,only_good)
%% User inputs:
% e.g. dataDir = 'Z:\Sherry\ephys_acquisition\SLV123\7turn_240731_133751';
% ksDir = 'Z:\Sherry\ephys_acquisition\SLV123\7turn_240731_133751\kilosort4';

%% to get spike times
nWF = 1e3;  % number of waveforms to get for each unit (1,000 spikes in)
spkOffset = 1.5e-3; % seconds before spike time to extract (1.5ms beforehand)
spkDur = 4e-3; % total duration(s) of waveform (4ms)
dt = 0.5/1e3; % time step for ccg binning (0.5ms)

%% extract intan raw data info 
h = readIntanInfo(fullfile(dataDir,'info.rhd')); %this is a very useful function that loads many intan infos. 
nCh = h.num_h.amplifier_channels;
fs = h.sample_rate; % the number of samples acquired in a second 

%% set up files
fDat = fullfile(dataDir,'amplifier.dat');
spkSamp = readNPY(fullfile(ksDir, 'spike_times.npy'));
sID = readNPY(fullfile(ksDir,'spike_clusters.npy'));
[unit_ID,cluster_labels] = getPhyClusterLabels(ksDir); % this loads the cluster quality rated from manual inspection in Phy
if only_good
    ind = strcmp(cluster_labels,'good'); % strcmp returns a list of logical values 
else
    ind = true(size(cluster_labels)); % all are 1 (true)
end
goodIDs = unit_ID(ind(:));
goodLabels = cluster_labels(ind(:));

%% Preparing variables
filenamestruct = dir(fDat);
nSamp = filenamestruct.bytes/(nCh*2);  % Number of samples per channel (int16 is 2 bytes each sample)
mmf = memmapfile(fDat, 'Format', {'int16', [nCh nSamp], 'x'}); % this creates a memory-mapped file organised in [nCh nSamp] matrix, and can be accessed by mmf.Data.x

numUnits = length(goodIDs);
spkDurSamples = spkDur*fs; % total samples to pre-allocate for waveform
maxTime = double(max(spkSamp))/fs; % the duration of recording (in this case 7800s)

% setting up matrices 
waveFormsMean = nan(spkDurSamples,nCh,numUnits); % create an empty matrix with the correct dimensions
meanRate = nan(numUnits,1); % the mean rate of each unit/cluster
medISI = nan(numUnits,1);
contam = nan(numUnits,1); % Why are we calculating contam here? HP
max_site = NaN(numUnits,1);
mxWF = nan(numUnits,spkDurSamples);
nSpikes = NaN(numUnits,1);

%% read in spikes for all units
for thisUnit=1:numUnits
    curUnitID = goodIDs(thisUnit);
    curSpikeTimes = double(spkSamp(sID==curUnitID))/fs;
    meanRate(thisUnit) = length(curSpikeTimes)./maxTime;
    medISI(thisUnit) = median(diff(curSpikeTimes));
    nSpikes(thisUnit) = length(curSpikeTimes);
    if isempty(curSpikeTimes); continue; end
    
    [K, Qi, Q00, Q01, Ri] = ccg(curSpikeTimes, curSpikeTimes, 500, dt); % compute autocorrelogram among spikes within a cluster
    contam(thisUnit) = min(Qi/(max(Q00, Q01)));
    % Q00:  measure of outer shoulder height, norm Poisson expectation
    % Q01: inner shoulder height, norm. by Poisson expectation
    % Qi: measure of inner refractory period height, different window
    % sizes (1,2,3 bins etc)
        
    % Exclude spikes whos waveform extends beyond recording
    curSpikeTimes((curSpikeTimes+spkDur*10)*fs > nSamp) = []; % spikes extending past recording
    curSpikeTimes((curSpikeTimes-spkOffset*10)*fs < 1) = []; % spikes starting before recording
    curUnitnSpikes = size(curSpikeTimes,1);
    
    spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
    spikeTimeKeeps = double(sort(spikeTimesRP(1:min([nWF curUnitnSpikes]))));
    
    mean_wave_filt = getSpikeWaveformNoFilter(mmf,spikeTimeKeeps, fs, spkOffset, spkDur); % filters removed by Sherry
    waveFormsMean(:,:,thisUnit) = mean_wave_filt;
    
    % Get channel with largest amplitude, take that as the waveform
    amps = max(mean_wave_filt)-min(mean_wave_filt);
    [max_val,max_site(thisUnit)] = max(amps); % Max site is the Intan channel number, so 1 = A-000
    mxWF(thisUnit,:) = waveFormsMean(:,max_site(thisUnit),thisUnit)';
    
    %         K((length(K)+1)/2) = 0;
    %         figure; subplot(2,2,1); plot(K); subplot(2,2,2); plot(Qi); title('Qi'); subplot(2,2,3); plot(Ri); title('Ri')
    %         subplot(2,2,4); plot(mxWF(thisUnit,:),'k')
    disp(['Completed ' int2str(thisUnit) ' units of ' int2str(numUnits) '.']);
    
end


