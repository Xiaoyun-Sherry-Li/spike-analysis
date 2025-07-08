function [spike_width, spike_asymmetry, firing_rate] = computeWaveformChar(wfs)
%% to find the troughs of waveforms
[trough,trough_idx] = min(wfs.mxWF, [], 2);% troughs happen at around time point 46 

% to find the two peaks before and after the trough 
peak_bef = zeros(size(wfs.mxWF,1),1);
peak_aft = zeros(size(wfs.mxWF,1),1);
peak_bef_time = zeros(size(wfs.mxWF,1),1);
peak_aft_time = zeros(size(wfs.mxWF,1),1);
for i = 1:size(wfs.mxWF,1)
    [peak_bef(i),peak_bef_time(i)] = max(wfs.mxWF(i,1:trough_idx(i)),[],2);
    [peak_aft(i),peak_aft_time(i)] = max(wfs.mxWF(i,trough_idx(i): trough_idx(i) + 30),[],2); % XL: note trough_idx(i) + 30 (1ms), instead of to the end, to avoid later irregularities in spike wf affecting determination of spike after trough peak  
end 

%%
spike_width =  peak_aft_time-1;

spike_asymmetry = zeros(size(wfs.mxWF,1),1);
for i = 1:size(wfs.mxWF,1)
    spike_asymmetry(i) = (peak_aft(i) - peak_bef(i))/(peak_aft(i) + peak_bef(i)); 
end
firing_rate = wfs.meanRate;