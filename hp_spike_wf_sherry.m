datatype = 'int16';
nCh = 64;
dataTypeNBytes = numel(typecast(cast(0, datatype), 'uint8'));
filenamestruct = dir('D:\ks_SLV123\amplifier.dat');
filenamestruct.bytes
nCh = 64;
nSamp = filenamestruct.bytes/(nCh*dataTypeNBytes);
mmf = memmapfile('D:\ks_SLV123\amplifier.dat', 'Format', {datatype, [nCh nSamp], 'x'});
spike_times = readNPY('Z:\Sherry\ephys_acquisition\SLV123\7turn_240731_133751\kilosort4\spike_times.npy');
spike_clusters = readNPY('Z:\Sherry\ephys_acquisition\SLV123\7turn_240731_133751\kilosort4\spike_clusters.npy');

idx209= find(spike_clusters == 209);
st_clust209 = double(spike_times(idx209));
Fs = 30000;
offset_t_final = 1000;
duration_t_final = 500;
getSpikeWaveform(mmf, st_clust209, Fs,offset_t_final, duration_t_final)