function [dig_raw, h, dig_time] = readIntanDig(filepath, tduration, tstart)
%READINTANDIG Reads Intan Technologies RHD2000 digital data file.
%% This is a modified version by Sherry (Feb 13, 2025) to account for the fact that in raw data two digital-ins were turned on but only one was used. 
%% the modification happened at line 63
%
%   [dig_raw, h] = READINTANDIG(filepath) reads the digital data from the
%   specified filepath. If filepath is not provided, a dialog will prompt
%   the user to select a folder. The function returns the raw digital data
%   in matrix form (n channels x n samples) and the header information.
%
%   [dig_raw, h] = READINTANDIG(filepath, tduration) reads the digital
%   data for a specified duration (in seconds) from the start of the file.
%
%   [dig_raw, h] = READINTANDIG(filepath, tduration, tstart) reads the
%   digital data for a specified duration (in seconds) starting from the
%   specified start time (in seconds).
%
%   [dig_raw, h, dig_time] = READINTANDIG(...) converts the output into timestamp format 
%
%   Inputs:
%       - filepath : String specifying the path to the data folder.
%                    If not provided, a dialog will prompt the user.
%       - tduration: Duration in seconds to read from the file.
%                    If not provided, the entire file is read.
%       - tstart   : Start time in seconds from the beginning of the file.
%                    If not provided, reading starts from the beginning.
%
%   Outputs:
%       - dig_raw : Matrix containing the raw digital data.
%                    Size is (n channels x n samples).
%       - h        : Struct containing header information from the info.rhd file.
%                    h.board_dig_in_channels.custom_channel_name contains
%                    channel names
%       - dig_time : Cell array, each element corresponds to one digital
%       input, and contains start (1st col) and stop (2nd col) 
%       timestamps for the digital input
%   Example:
%        [dig_raw, h, dig_time] = readIntanDig;
%       Prompts for a raw folder and provides all outputs for the entire file.
%       [dig_raw, h] = readIntanDig('D:\data\HC15_231007\raw_231007_125517', 10, 5);
%       This reads 10 seconds of data starting from the 5th second.
%
%   Note:
%       This function is based on read_Intan_RHD2000_file and has been
%       edited by Hannah Payne, Aronov lab.
%
%   See also: READINTANINFO, READINTANAMP

if ~exist('filepath','var') || isempty(filepath)
    [filepath] = ...
        uigetdir('Select a folder with single file per data type data');
    if (filepath == 0);  return; end
end

% Open the info file
info_filepath = fullfile(filepath,'info.rhd');
fprintf('Opening info file %s\n',info_filepath)
h = readIntanInfo(info_filepath);
% ndig_in = h.num_h.board_dig_in_channels;  % Sherry commented it out
% because in the acquisition file she had to digital-in channels on, but
% she only used the first one. Without specifying the first one throws an
% error later on. 
ndig_in = 1; % sherry added this 
dig_raw = [];

if ndig_in
    % Open the digital data
    dig_filepath = fullfile(filepath,'digitalin.dat');
    fid = fopen(dig_filepath, 'r');
    
    % Get the total number of time point samples
    fileinfo = dir(dig_filepath);
    num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
    
    % Change start time if needed
    if exist('tstart','var') && ~isempty(tstart)
        ind_start = round(tstart*h.sample_rate);
        fseek(fid, ind_start*2, 'bof');
        num_samples = num_samples - ind_start;
    else
        tstart = 0;
    end
    
    % Change the duration to read if needed
    if exist('tduration','var') && ~isempty(tduration)
        num_samples = round(tduration*h.sample_rate);
    end
    
    digital_word = fread(fid, num_samples, 'uint16');
    fclose(fid);
    
    % Individual digital inputs can be isolated using the bitand function in MATLAB:
    % digital_input_ch = (bitand(digital_word, 2^ch) > 0); % ch has a value of 0-15 here
    dig_raw = false(ndig_in, length(digital_word));
    chanlabel = cell(ndig_in,1);

    for ii = 1:ndig_in
        ch = h.board_dig_in_channels(ii).native_order;
        chanlabel{ii}  = h.board_dig_in_channels(ii).custom_channel_name;
        dig_raw(ii,:) = (bitand(digital_word, 2^ch) > 0);
    end
end


if nargout>2
    samplerate = h.sample_rate;
    nt = size(dig_raw,2);
    t_dig = (0:nt-1)/samplerate + tstart; % (s)
    dig_time = cell(ndig_in,1);
    for ii = 1:ndig_in
        iData_rise = t_dig(diff(dig_raw(ii,:))>0);
        iData_fall = t_dig(diff(dig_raw(ii,:))<0);
        iData_rise = iData_rise(:);
        iData_fall = iData_fall(:);
        
        % First column is rise times, second column is fall times
        if iData_rise(1)>iData_fall(1)
            iData_rise = [NaN; iData_rise];
        end
        if iData_fall(end)<iData_rise(end)
            iData_fall = [iData_fall; NaN];
        end
        dig_time{ii} = [iData_rise iData_fall];
    end
end
fprintf('Read digital data from channel %s\n',chanlabel{:})
