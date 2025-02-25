function [map, info, bin_centers, pos_time_smooth] = placeCellAnalysis(...
    xx, yy, fps, ispk, nshuffle, lims, nbins, nhamming, min_occup, plot_on)
% PLACECELLANALYSIS Create place maps and calculate spatial info
%
% map = placeCellAnalysis(xx, yy, fps, ispk, nshuffle, lims, nbins, nhamming, min_occup)
% calculates the place map. Inputs:
%   XX and YY: animal's position *(e.g. in mm, but so long as you are consistent throughout any units ok)
%   FPS: frames per second of the behavior
%   ISPK: spike indices into the behavior ISPK
%   NSHUFFLE: number of random shuffles to conduct (specify 0 for no shuffling)
%   LIMS: two-number vector limits of the square arena (*[0 609.6] (mm) in the 2-foot arena)
%   NBINS: number of bins in which to divide the behavior *(40 bins)
%   NHAMMING: size of hamming filter *(13)
%   MIN_OCCUP: minimum occupancy time in each bin *(0.1 s)
%
% map = placeCellAnalysis(..., plot_on) plots the results
%
% [map, info] = placeCellAnalysis(...) also returns spatial information
% according to Skaggs et al. 1993, in bits/spike.
%
% [map, info, bin_centers] = placeCellAnalysis(...) also returns the bin
% centers corresponding to the place map
%
% [map, info, bin_centers, pos_time_smooth] = placeCellAnalysis(...) also
% returns the smoothed occupancy map used to calculate the place map
%
% *parameters used in Payne et al. 2021
%
% Dependencies: nanconv
%
% Hannah Payne, Aronov Lab 2018


if ~exist('plot_on','var')
    plot_on = 0;
end

% Preallocate info vector
info = zeros(1,max(1,nshuffle));

% Get bin edges
bin_edges = linspace(lims(1), lims(2), nbins+1);
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;

% Get hamming filter
w = hamming(nhamming)*hamming(nhamming)';
w = w/sum(w(:));

% Calculate position counts (doesn't change with shuffling)
% occupancy map
pos_count = histcounts2(xx, yy, bin_edges, bin_edges)'; % Transpose because it places the x categories in rows not columns
pos_time = pos_count/fps;
pos_time_smooth = nanconv(pos_time, w,'nanout', 'edge'); % no difference from conv if there aren't any nans, but important when making gaze maps

map = NaN(nbins, nbins, nshuffle);
for i_shuffle = 1:max(1,nshuffle)
    if nshuffle > 0
        r = ceil(rand*length(xx));
        ispk = mod(ispk + r - 1, length(xx)) + 1;        % faster than alternatice
    end
    
    % Get xy locations during spike times
    xx_sp = xx(ispk);
    yy_sp = yy(ispk);
    
    % Get heat map
    sp_count = histcounts2(xx_sp,yy_sp, bin_edges, bin_edges)';
    
    % Filter while dealing with NaNs (file exchange function nanconv)
    sp_count_smooth = nanconv(sp_count, w,'nanout', 'edge');
    curr_map = sp_count_smooth./pos_time_smooth;
    
    % Require a minimum time spent in this region
    curr_map(pos_time_smooth<=min_occup) = NaN;
    pos_time_smooth(pos_time_smooth<=min_occup) = NaN;
    
    % Spatial information per spike (note: this is bits/spike, not bits per
    % second like in Skaggs et al. 1993).
    f = ~isnan(curr_map);
    h = sp_count_smooth(f);
    t = pos_time_smooth(f);
    lambdai = h./t;             % should equal map(f)
    p_i = t/sum(t);
    lambda = sum(p_i.*lambdai); % mean rate
    % lambda = sum(h)/sum(t); % same as above
    f = lambdai>0;
    info(i_shuffle) = sum(p_i(f).*(lambdai(f)/lambda).*log2(lambdai(f)/lambda));
    map(:,:,i_shuffle) = curr_map;
    
end


if plot_on
    
    figure;
    ah = subplot(121);
    ah(2) = subplot(122);
    
    % Plot behavioral trajectory
    h = plot(ah(1), xx,yy,'color',[0 0 0],'LineWidth',.25); hold(ah(1),'on');
    
    % Plot red dots for spikes
    if length(ispk)<4000;    msize = 3;  else;      msize= 1;    end
    h(2) = plot(ah(1), xx(ispk),yy(ispk),'or','MarkerSize',msize,'MarkerFaceColor','r','MarkerEdgeColor','none');
    
    % Plot heat map
    h(3) = imagesc(ah(2), bin_centers([1 end]),bin_centers([1 end]),curr_map,'AlphaData',~isnan(curr_map));
    clims = [0 max(1, prctile(curr_map(:),99))];
    set(ah(2),'CLim', clims)
    set(ah,'YDir','normal')
    set(ah,'XTick',[],'YTick',[])
    axis(ah,'image');
    set(ah,'XLim',bin_edges([1 end]),'YLim',bin_edges([1 end]));
    drawnow

end


