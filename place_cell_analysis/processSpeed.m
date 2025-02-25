function B = processSpeed(B, tmaxgap, tsmooth_linear_speed, speed_thresh, max_duration, ploton)
% PROCESSSPEED - Process and filter head speed data
%.
0
% Usage:
%   B = processSpeed(B, tmaxgap, tsmooth_linear_speed, speed_thresh, max_duration, ploton)
%
% Inputs:
%   B                     - Structure containing tracking data with fields:
%                           B.head_xpos, B.head_ypos, B.head_zpos, and B.tt (timestamps).
%   tmaxgap               - Maximum gap (in seconds) allowed for interpolation.
%   tsmooth_linear_speed  - Time window (in seconds) for smoothing speed.
%   speed_thresh          - Speed threshold (mm/s) for masking low-speed segments.
%   max_duration          - Maximum duration (seconds) of stationary periods to exclude.
%   ploton (optional)     - If true, plots speed before and after thresholding.
%
% Outputs:
%   B                     - Updated structure with processed speed and filtered positions.
%
% Description:
%   - Interpolates missing head position data within the specified gap threshold.
%   - Computes smoothed linear speed based on head movements in the XY plane.
%   - Identifies and masks low-speed periods based on threshold and duration.
%   - Updates head position and speed fields, setting masked values to NaN.
%   - Optionally plots the speed before and after masking.
%
% Notes:
%   - Speed is computed in mm/s and smoothed using a moving average filter.
%   - Masking is applied to head_xpos, head_ypos, head_zpos, and other relevant fields.
%
% Example:
%   B = processSpeed(B, 0.2, 1, 50, 5, true);
% 
% Parameters from Payne et al. 2021:
%   speed_thresh = 50;        % (mm/s)
%   max_duration = 5;         % (s)
%   tsmooth_linear_speed = 1; % (s)
%   tmaxgap = 0.2;            % (s) 

% Get the frame rate
fps = double(1/mean(diff(B.tt)));

% Fill in small gaps in the behavior data
% nmaxgap = round(tmaxgap*fps);
% B.head_xpos = interp1gap(B.head_xpos, nmaxgap);
% B.head_ypos = interp1gap(B.head_ypos, nmaxgap);
% B.head_zpos = interp1gap(B.head_zpos, nmaxgap); % optional, delete if unused

% Get the linear speed of the head (in mm/s) projected on xy plane
% nsmooth_linear_speed = round(tsmooth_linear_speed*fps/2)*2+1;
% smooth_xx = smooth(B.head_xpos, nsmooth_linear_speed , 'moving');
% smooth_yy = smooth(B.head_ypos, nsmooth_linear_speed , 'moving');
% smooth_zz = smooth(B.head_zpos, nsmooth_linear_speed , 'moving'); % sherry added the z axis
% B.speed = [sqrt(diff(smooth_xx).^2+diff(smooth_yy).^2 + diff(smooth_zz).^2)*fps; NaN];

% sherry removed smoothing (already done in python) 
B.speed = [sqrt(diff(B.head_xpos).^2+diff(B.head_ypos).^2 + diff(B.head_zpos).^2)*fps; NaN];

% Create mask for speed thresholding
if speed_thresh
    if max_duration
        code_area = regionprops(B.speed<speed_thresh,'area');
        code_label = bwlabel(B.speed<speed_thresh);
        stationary_lengths = [code_area.Area]/fps;
        reject_segs = find(stationary_lengths > max_duration);
        B.speed_mask = false(size(B.tt));
        for jj = 1:length(reject_segs)
            B.speed_mask(code_label==reject_segs(jj)) = 1;
        end
    else
        B.speed_mask = B.speed < speed_thresh;
    end
    
    % Plot speed thresholding before and after
    if exist('ploton','var') && ploton
        figure; plot(B.tt,B.speed,'k'); hold on
        temp = B.speed; temp(B.speed_mask) = NaN;
        plot(B.tt, temp,'r');   xlabel('Time (s)'); ylabel('Speed (mm/s)')
    end
    fprintf('%f excluded by speed thresholding\n', mean(B.speed_mask))
    
end

% Apply the speed thresholding to the output
B.speed_raw = B.speed; % Keep a copy to analyze speed tuning
B.head_xpos(B.speed_mask) = NaN;
B.head_ypos(B.speed_mask) = NaN;
B.head_zpos(B.speed_mask) = NaN;
B.speed(B.speed_mask) = NaN;

