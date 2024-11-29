function [PD_onset_tr] = get_PD_trial(Z_ML, Z_Digital, Z_Analog)

%   get_PD_trial.m gets Z files and passes the Photodiode onsets in secs
%
%   Inputs:
%       Z Files
%
%   Outputs:
%       PD_onset_tr: n*3 matrix, n = number of trials, the number is in sec
%       from beginning of the file
%       
   
    st_starts = Z_Digital.data(2, find(Z_Digital.data(1, :) == 9));
    st_ends = Z_Digital.data(2, find(Z_Digital.data(1, :) == 18));
                   
    PD_Threshold_mV = 2000;
    PD = Z_Analog.data(4, :);
    PD_timestamps = (0:length(PD)-1) ./ Z_Analog.info.SampleRate;
    PD_onset_tr = nan(length(Z_ML.condition), 3);
    
    for j = 1:length(Z_ML.condition)
        PD_tr = PD(find(PD_timestamps >  st_starts(j) & PD_timestamps <  st_ends(j)));
    
        PD1 = find(diff(PD_tr) > PD_Threshold_mV, 1); % in sample
        if(~isempty(PD1)); PD_onset_tr(j, 1) = st_starts(j) + PD1/Z_Analog.info.SampleRate; end
        PD2 = find(diff(PD_tr) > PD_Threshold_mV, 2); % in sample
        if(~isempty(PD2) && length(PD2)>= 2 && PD2(2) > PD2(1) + 10); PD_onset_tr(j, 2) = st_starts(j) + PD2(2)/Z_Analog.info.SampleRate; end
        PD3 = find(diff(PD_tr) > PD_Threshold_mV, 3); % in sample
        if(~isempty(PD3) && length(PD3)>= 3 && PD3(3) > PD3(2) + 10 && PD3(2) > PD3(1) + 10); PD_onset_tr(j, 3) = st_starts(j) + PD3(3)/Z_Analog.info.SampleRate; end
        
    end

end