dataDir = 'Z:\Sherry\ephys_acquisition\SLV123\7turn_240731_133751';
ksDir = "Z:\Sherry\ephys_acquisition\SLV123\7turn_240731_133751\kilosort4\";
wfs = getSessionWaveforms(dataDir, ksDir);

wfNoFilter = extractWaveforms(dataDir, ksDir,true);


%% to find the troughs of waveforms
[trough,trough_idx] = min(wfs.mxWF, [], 2);% troughs happen at around time point 46 
% to find the two peaks before and after the trough 
peak_bef = zeros(size(wfs.mxWF,1),1);
peak_aft = zeros(size(wfs.mxWF,1),1);
peak_bef_time = zeros(size(wfs.mxWF,1),1);
peak_aft_time = zeros(size(wfs.mxWF,1),1);
for i = 1:size(wfs.mxWF,1)
    [peak_bef(i),peak_bef_time(i)] = max(wfs.mxWF(i,1:trough_idx(i)),[],2);
    [peak_aft(i),peak_aft_time(i)] = max(wfs.mxWF(i,trough_idx(i):120),[],2);
end 

spike_width =  peak_aft_time-1

spike_asymmetry = zeros(size(wfs.mxWF,1),1);
for i = 1:size(wfs.mxWF,1)
    spike_asymmetry(i) = (peak_aft(i) - peak_bef(i))/(peak_aft(i) + peak_bef(i));
end
firing_rate = wfs.meanRate;
wf_chars = [spike_asymmetry,spike_width,firing_rate];

% Perform k-means clustering
k = 3
[idx, C] = kmeans(wf_chars, k);

% Visualize the clustering
figure;
scatter3(wf_chars(:,1), wf_chars(:,2), wf_chars(:,3), 50, idx, 'filled');
set(gca, 'ZScale', 'log')
hold on;
scatter3(C(:,1), C(:,2), C(:,3), 100, 'kx', 'LineWidth', 2);
title('3D Cluster Visualization');
xlabel('spike_asymmetry');
ylabel('spike_width');
zlabel('firing_rate');
grid on;

%% to further examine the extremely high firing rate interneurons (cluster 3)
clust3_idx = find(idx == 3);
meanrate = zeros(4,1);
for i = 1:4
    plot(wfs.mxWF(clust3_idx(i),:))
    grid on
    meanrate(i) = wfs.meanRate(clust3_idx(i));
    %text(5, 100+100i, ['FR =' num2str(meanrate)]);
    %text(5*i, 300, ['Cell ID = ' num2str(clust3_idx(1))]);
    hold on 
end 



