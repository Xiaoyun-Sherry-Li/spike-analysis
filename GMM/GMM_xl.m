clear all
close all 
clc
%%
ephyspath = 'Z:\Sherry\acquisition\LIM124_082125\LIM124_8turn_250821_100934_v1';
wfs = getSessionWaveforms(ephyspath, fullfile(ephyspath,'kilosort4'),1);

%%
save(fullfile(path,'wfs_LIM124_082125_v1.mat'),"wfs");
% [spike_width, spike_asymmetry, firing_rate] = computeWaveformChar(int.wfs);

%%
path = 'Z:\Sherry\ephys_acquisition\WFstruct';
wfs_list = dir(fullfile(path,'wfs_LVN4_041125*'));

%%
[gm, idx,labels] = GMM_make(wfs_list, 1, 2);

%% 