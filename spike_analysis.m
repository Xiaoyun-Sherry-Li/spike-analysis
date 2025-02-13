% load data properties
KS_dir = 'Z:\Sherry\ephys_acquisition\SLV123\7turn_240731_133751\kilosort4';
sp = loadKSdir(KS_dir); 
spike_times = readNPY('Z:\Sherry\ephys_acquisition\SLV123\7turn_240731_133751\kilosort4\spike_times.npy');
% are the spike times align to the peak of a spike? 
spike_clusters = readNPY('Z:\Sherry\ephys_acquisition\SLV123\7turn_240731_133751\kilosort4\spike_clusters.npy');
cluster_group = readtable('Z:\Sherry\ephys_acquisition\SLV123\7turn_240731_133751\kilosort4\cluster_group.tsv','FileType', 'text', 'Delimiter', '\t');

%fseek, fread, fclose, 
% first, try to locate the time points when spikes happen

% to determine the size of raw data file in uint8 (one bit):
datatype = 'int16';
dataTypeNBytes = numel(typecast(cast(0, datatype), 'uint8'));
filenamestruct = dir('D:\ks_SLV123\amplifier.dat');
filenamestruct.bytes
nCh = 64;
nSamp = filenamestruct.bytes/(nCh*dataTypeNBytes);

% to extract at the first spike position (t = 2 * nCh = 64)
mmf = memmapfile('D:\ks_SLV123\amplifier.dat', 'Format', {datatype, [nCh nSamp], 'x'});

% how to pick the max chan? probably the median value of each channel? 
% amplifier_data = 0.195*mmf.Data.x(3,3:10);

fid = fopen('D:\ks_SLV123\amplifier.dat');

% to 1) locate the spike times and its corresponding clusters, extract  
spike_time = 326;
spike_cluster = spike_clusters(spike_time);
cluster_quality = cluster_group.group(cluster_group.cluster_id == spike_cluster);
% around each spike times, which should be where the peak of the spike is,
% to extract +/-1ms snapshoots, which is 2*30 = 60 timepoints
spike_loc_2ms = (spike_times(spike_time)-30)*64; %temporaily changed to 1000 instead of 60
% to convert 120 timepoints duration into number of bytes: the data type is
% int16, so each data point is 2 bytes, and each time point has 64
% channels, and there are 120 time points: 64 * 120 = 7680
four_ms = 30 * 2;%temporaily changed to 90 instead of 4
num_points_4ms = nCh * four_ms;
fseek(fid,spike_loc_2ms,'bof');
st8 = fread(fid,num_points_4ms,'*int16');
chan_time = reshape(st8,64,60);
chan_time_amp = chan_time * 0.195;

% to map channel positions using chanMap %% sorting order based on chanMap is very tricky
importfile("D:\GitHub\spike_sort\spikesort-hp\src\kilosort_config_files\H5.mat");
[~,idx] = sort(chanMap);
sorted_chan_time = chan_time_amp(idx,:);
median_chans = median(st8,2);
sorted_st8_med = round(median_chans(idx));

% in this case, chan 31 in sorted_st3 is largest

y= -300:300;
hold on;
for i = 1:64
    plot(sorted_chan_time(i,:));
end

y= -300:300;
plot(sorted_chan_time(62,:));
hold on;
for i = 1:5
    plot(sorted_chan_time(i,:));
end
hold on;

% chan 20 and 53 should be broken channels, according to intan, should have
% strong waves
% but from median values, chan 30 and 31 seems to be problematic (at st 10,
% sc 53), chan 11 and 16 (at st 8, sc 18), chan 30 and 31 again at st 15,
% sc 209, chan 30 and 31 again at st 40, sc 132. chan 11 and 16 at st 70,
% chan 30 and 31 at st 5000, sc 3
% sc 170.

idx209= find(spike_clusters == 209);
st_clust209 = spike_times(idx209);
st_clust209(1:100);

fid = fopen('D:\ks_SLV123\amplifier.dat');

% to 1) locate the spike times and its corresponding clusters, extract  
spike_time = st_clust209;
%spike_cluster = spike_clusters(spike_time(1));
cluster_quality = cluster_group.group(cluster_group.cluster_id == 209);


% around each spike times, which should be where the peak of the spike is,
% to extract +/-1ms snapshoots, which is 2*30 = 60 timepoints
spike_loc_2ms = (spike_times(spike_time(1))-30)*64; %temporaily changed to 1000 instead of 60
% to convert 120 timepoints duration into number of bytes: the data type is
% int16, so each data point is 2 bytes, and each time point has 64
% channels, and there are 120 time points: 64 * 120 = 7680
four_ms = 30 * 2;%temporaily changed to 90 instead of 4
num_points_4ms = nCh * four_ms;
fseek(fid,spike_loc_2ms,'bof');
st8 = fread(fid,num_points_4ms,'*int16');
chan_time = reshape(st8,64,60);
chan_time_amp = chan_time * 0.195;

% to map channel positions using chanMap %% sorting order based on chanMap is very tricky
importfile("D:\GitHub\spike_sort\spikesort-hp\src\kilosort_config_files\H5.mat");
[~,idx] = sort(chanMap);
sorted_chan_time = chan_time_amp(idx,:);
median_chans = median(st8,2);
sorted_st8_med = round(median_chans(idx));

% in this case, chan 31 in sorted_st3 is largest

y= -300:300;
hold on;
for i = 1:64
    plot(sorted_chan_time(i,:));
end

y= -300:300;
plot(sorted_chan_time(62,:));
hold on;
for i = 1:5
    plot(sorted_chan_time(i,:));
end


fclose(fid);
