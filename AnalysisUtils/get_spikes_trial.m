function [spk_mat, spdn_mat] = get_spikes_trial(Z_Analog, spikes, tr, varargin)

%   get_spikes_trial.m gets a spike data and passes the aligned continous data
%
%   Inputs:
%       ContData:   The continous data that you wanna split
%       spikes:     The spike channel, an array of time points from start
%                   of the recording, this is one spike array so choose the unit and
%                   channel before this function
%       tr:         list of trials that you want the splited data
%       *PDnumber:  which photodiode onset you want to use for alignment
%                   (default: 1, this is an optional input when you have more than 1
%                   photodiode onset on a trial)
%       *predur:    duration in ms before the Photodiode onset
%       *postdur:   duration in ms after the Photodiode onset
%   
%   Outputs:
%       spk_mat:    Matrix of spikes in size N*m which N is the number of trials
%                   that you want and m is the length you determined using pre_dur and
%                   post_dur, m = pre_dur + post_dur + 1
%       spdn_mat:   Matrix of spdn in size N*m which N is the number of trials
%                   that you want and m is the length you determined using pre_dur and
%                   post_dur, m = pre_dur + post_dur + 1
%
%   Notes: takes ~ 100 ms for 1000 trials
%
%   Sample:         
%           cond_want = 5;          %5: ETP
%           trs_right = find(Z_ML.condition == cond_want & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Right'));
%           spk = Z_Spikes.data(ch).neural_timeStamps(1, find(Z_Spikes.data(ch).neural_timeStamps(2, :) == Z_Spikes.data(ch).Units(unitID)));
%           [SPIKES_RIGHT, SPDEN_RIGHT] = get_spikes_trial(Z_Analog, spk, trs_right);

    
    p = inputParser;
    p.addOptional('PDnumber', 1);
    p.addOptional('predur', 500);
    p.addOptional('postdur', 500);
    p.parse(varargin{:});
    
    pre_dur = p.Results.predur;
    post_dur = p.Results.postdur;

    % ____________________________ SPDN part ______________________________
    % The spden function will be created by the EPSP function. Here is how it
    % is defined
    growth = 1;
    decay = 20;
    A = [];
    for i = 0:1:100  % A 100 ms function
        A(end+1,1) = i;
        A(end,2) = (1 - exp(i*-1/growth)) .* (exp(i*-1/decay));
    end
    B = cumsum(A,1);    %Sums along columns to extract area
    C = B(size(B,1),size(B,2)); % Value in last row and column of B is the area of the curve
    KERNEL = [A(:,1) A(:,2)/C];
    % _____________________________________________________________________
               
    time_range = -pre_dur:post_dur;
    spk_mat = zeros(length(tr), length(time_range));

    switch p.Results.PDnumber
        case 1
            PD = Z_Analog.PD_sec(:, 1);
        case 2
            PD = Z_Analog.PD_sec(:, 2);
        case 3
            PD = Z_Analog.PD_sec(:, 3);
    end

    CC = zeros(length(tr), length(time_range)); % Row of zeros for timespan of trial

    for j = 1:length(tr)

        if(isnan(PD(tr(j))))
            warning('Not enough PD crossings in trial %d.', tr(j));
            continue;
        end
        
        spike_times = (spikes(spikes > PD(tr(j)) - pre_dur/Z_Analog.info.SampleRate & spikes < PD(tr(j)) + post_dur/Z_Analog.info.SampleRate) - PD(tr(j)))*1000;
        spk_mat(j, 1:length(spike_times)) = spike_times;
        CC(j, round(spike_times+pre_dur)+1) = 1;
    end

    spdn_mat = conv2(CC, KERNEL(:,2)') * 1000;
    spdn_mat = spdn_mat(:, 1:length(time_range));

end