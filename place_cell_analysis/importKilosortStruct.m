function S = importKilosortStruct(ksDir, fs, option_only_good, session_label)
% IMPORTKILOSORTSTRUCT import sorted spike times and waveforms and convert to dat
% structure
%
% S = IMPORTKILOSORTSTRUCT(filepath, fs) imports all available channels of sorted data
% 	from the folder filepath, recorded at samplerate fs (usually 30000)
% option_only_good: Set 1 to only load "good" cells (post manual phy sorting, if available)
% session_label: a prefix added to each cell name (before K000 assigned unit number etc)
%
% Example:
% S = importKilosortStruct('Z:\Hannah\ephys\project2\HC05_220825\kilosort2_output', 3e4, 1,'HC05_220825')

% read in spike times and get IDs and kilosort labels
tm = readNPY(fullfile(ksDir,'spike_times.npy'));
tm = double(tm)/fs;
sID = readNPY(fullfile(ksDir,'spike_clusters.npy'));
% [unit_ID, goodLabels] = get_phy_cluster_labels(ksDir);
[unit_ID, goodLabels] = getPhyClusterLabels(ksDir);

% Load these from waveformStruct.mat!
wvStruct = getfield(load(fullfile(ksDir, 'waveformStruct.mat')),'wvStruct');
if ~all(unit_ID(:) == wvStruct.goodIDs(:)); warning('Check spike sorting'); end


% Only look at "good" clusters
if option_only_good
    inds_keep = find(strcmp(goodLabels,'good') | goodLabels==2);
else % Discard "noise" sorted clusters
    inds_keep = find(~strcmp(goodLabels,'noise') | goodLabels>0);
end
inds_keep = inds_keep(:)';

% Loop through units.
S = struct;
ii_out = 1;
for ii = inds_keep(:)'
    
    % Cell ID
    cID = unit_ID(ii);
    
    % Get spike times from this unit
    tt = tm(sID==cID);
    
    % Store spikes
    if exist('session_label','var')&& ~isempty(session_label)
        chanlabel =  sprintf('%s_K%03i',session_label, cID);
    else
        chanlabel =  sprintf('K%03i', cID); 
    end
    tstart = min(tt);
    tend = max(tt);
    
    % Include waveform and store as dat structure
    maxsite = wvStruct.max_site(ii); % intan channel index of largest waveform. So 1 = A-000.
    waveform = wvStruct.mxWF(ii,:); % mxWF stores largest waveform!
    info.waveunit = 'uV';
    info.wavefreq = fs;
    info.prethresh = fs*wvStruct.spkOffset;
    info.goodLabel = wvStruct.goodLabels{ii};
    
    % Load GMM results if available
    if isfield(wvStruct,'typeLabels')
        info.typeLabel = wvStruct.typeLabels{ii};
    end
    
    % Store this cell
    S(ii_out).data = tt;
    S(ii_out).chanlabel = chanlabel;
    S(ii_out).maxsite = maxsite;
    S(ii_out).samplerate = 'event';
    S(ii_out).tstart = tstart;
    S(ii_out).tend = tend;
    S(ii_out).units = 's';
    
    S(ii_out).cellID = cID;
    S(ii_out).waveform = waveform;
    S(ii_out).info = info;
    
    ii_out = ii_out+1;
end


