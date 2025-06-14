clc
clear all 
close all
addpath 'C:\Users\xl313\OneDrive\Documents\GitHub\spike_sort_Xl\spikesort-hp\src\kilosort_utils'
addpath 'C:\Users\xl313\OneDrive\Documents\GitHub\spike_sort_Xl\Kilosort-2.0\postProcess'
addpath 'C:\Users\xl313\OneDrive\Documents\GitHub\spike-analysis\extractWF_xl'
addpath 'C:\Users\xl313\OneDrive\Documents\GitHub\spike-analysis\npy-matlab-master'

dataDir = 'C:\Users\xl313\OneDrive\Desktop\3TurnsLVN4_250411_111802';
ksDir = fullfile(dataDir,'kilosort4');
wfs = getSessionWaveforms(dataDir, ksDir);

% wfNoFilter = extractWaveforms(dataDir, ksDir,true); temporaily commented
% out by Sherry 

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
%%
spike_width =  peak_aft_time-1;

spike_asymmetry = zeros(size(wfs.mxWF,1),1);
for i = 1:size(wfs.mxWF,1)
    spike_asymmetry(i) = (peak_aft(i) - peak_bef(i))/(peak_aft(i) + peak_bef(i));
end
firing_rate = wfs.meanRate;
wf_chars = [spike_asymmetry,spike_width,log10(firing_rate)];

% Perform k-means clustering. Remember to 1) log all firing rates before k means clustering 2) z score all data
k = 2;
[idx, C] = kmeans(zscore(wf_chars), k);
% Define colors for clusters
colors = [1 0 0; 0 0 1;0 1 0]; % Blue and red

% Map cluster indices to colors
cluster_colors = colors(idx, :);

% Visualize the clustering
figure;
scatter3(wf_chars(:,1), wf_chars(:,2), 10.^wf_chars(:,3), 25, cluster_colors);

set(gca, 'ZScale', 'log')
hold on;
% scatter3(C(:,1), C(:,2), C(:,3), 100, 'kx', 'LineWidth', 2);
title('K-means clustering');
xlabel('Spike asymmetry');
ylabel('Spike width');
zlabel('Firing rate');
grid on;

%% [optional] to further examine the extremely high firing rate interneurons (cluster 3)
t  = -39:80;
t = t/30;

clust3_idx = find(idx == 2);
meanrate = zeros(4,1);
for i = 1:length(clust3_idx)-1
    plot(t,wfs.mxWF(clust3_idx(i),:),'blue')
    grid on
    meanrate(i) = wfs.meanRate(clust3_idx(i));
    %text(5, 100+100i, ['FR =' num2str(meanrate)]);
    %text(5*i, 300, ['Cell ID = ' num2str(clust3_idx(1))]);
    hold on 
end 

hold on  

%% plot each neuron's mean wfs 
t  = -39:80;
t = t/30;

meanrate = zeros(4,1);
for i = 1:length(wfs.meanRate)
    figure;
    hold on
    plot(t,wfs.mxWF(i,:),'red')
    grid on
    meanrate(i) = wfs.meanRate(i);
    title(sprintf('Neuron ID: %d, FR = %.2f Hz', i, meanrate(i)));
    filename = fullfile('Z:\Sherry\acquisition\LVN4_041125\wfs_goodunits',sprintf('neuron%d.png', i));
    saveas(gcf,filename);
end 

hold off
