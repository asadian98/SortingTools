function [cont_tr] = get_continuous_trial(Z_Analog, ContData, tr, varargin)

%   get_continuous_trial.m gets a continous data and passes the aligned continous data
%
%   Inputs:
%       ContData:   The continous data (LFP, EMG, Eyedata) that you wanna
%       split (n * m) n: number of channels 
%       pre_dur:    duration in ms before the Photodiode onset
%       post_dur:   duration in ms after the Photodiode onset
%       ch:         channel number for the input signal
%       tr:         list of trials that you want the splited data
%       *PDnumber:  which photodiode onset you want to use for alignment
%       (default: 1, this is an optional input when you have more than 1
%       photodiode onset on a trial)
%
%   Outputs:
%       cont_tr:    Matrix in size K*N*m which K is the number of channels, N is the number of trials
%       that you want and m is the length you determined using pre_dur and
%       post_dur, m = pre_dur + post_dur + 1
%
%   Notes:
%       Default fs = 1k, pass fs as input for EMG data if you have 30k,
%       'fs', 30000
%       Default predur = 500 ms, postdur = 500 ms, pass them as input if
%       you want different values
%       'predur', 300
%       'postdur', 1000
%
%   takes ~60 for 1000 trials
%   It will give you a zero array if your trial doesn't have PD onset
%
%   TODO: check for EMG
% 

    p = inputParser;
    p.addOptional('PDnumber', 1);
    p.addOptional('fs', 1000);
    p.addOptional('predur', 500);
    p.addOptional('postdur', 500);
    p.parse(varargin{:});
    
    fs = p.Results.fs;
    pre_dur = p.Results.predur;
    post_dur = p.Results.postdur;
    
    cont_tr = zeros(size(ContData, 1), length(tr), pre_dur+post_dur+1);
    switch p.Results.PDnumber
        case 1
            PD = Z_Analog.PD_sec(:, 1);
        case 2
            PD = Z_Analog.PD_sec(:, 2);
        case 3
            PD = Z_Analog.PD_sec(:, 3);
    end

    for j = 1:length(tr)

        if(isnan(PD(tr(j))))
            warning('Not enough PD crossings in trial %d.', tr(j));
            continue;
        end
    
        cont_tr(:, j, :) = ContData(:, (round(PD(tr(j))*fs) - pre_dur*(fs/Z_Analog.info.SampleRate)):(round(PD(tr(j))*fs) + post_dur*(fs/Z_Analog.info.SampleRate)));
    end

end