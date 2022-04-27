function [locs, stim_interval] = findStims(trace)

% Description:
% Searches for high-amplitude stimulation artifacts.
%
% INPUT:
% trace [vector] - filtered voltage trace to be analysed
%
% OUTPUT:
% locs [vector]          - recording frames where artifacts were detected
% stim_interval [scalar] - estimated interval between subsequent stimuli

% Author:
%   Jeremy Chabros, University of Cambridge, 2021
%   email: jjc80@cam.ac.uk
%   github.com/jeremi-chabros

%%

recording_duration = length(trace);
fs = 25000;
thr = std(trace);

stim_interval = [];

[~,locs] = findpeaks(trace,...
    'MinPeakDistance',0.2*fs,...
    'Threshold',2*thr,...
    'MinPeakHeight',4*thr);

if length(locs) > 5
    t_vec = locs(1:end-1) - locs(2:end);
    stim_interval = -mode(t_vec);
    stim_fr = 1/round(stim_interval/fs);
    first_stim_time = locs(2); % does not have to be 1st, let's fill from both sides
    n_stims_possible = ceil(recording_duration/stim_interval);

    if n_stims_possible ~= length(locs)

        fill_val = first_stim_time;
        stims_filled = [];
        counter = 1;

        while fill_val > 0
            if counter == 1
                fill_val = first_stim_time - stim_interval;
                if fill_val > 0
                    stims_filled(counter) = fill_val;
                end
            else
                fill_val = stims_filled(counter-1) - stim_interval;
                if fill_val > 0
                    stims_filled(counter) = fill_val;
                end
            end
            counter = counter+1;
        end

        counter = length(stims_filled)+1;
        fill_val = max(stims_filled);

        while fill_val < recording_duration

            fill_val = stims_filled(counter-1) + stim_interval;
            if fill_val <= recording_duration
                stims_filled(counter) = fill_val;
            end
            counter = counter+1;
        end

        locs = stims_filled;
    end
end
if ~stim_interval
    stim_interval = 0;
end
end