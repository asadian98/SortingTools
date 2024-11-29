% This script demonstrates how to use the AnalysisUtils functions. 
% These versatile functions allow you to work with Z files and visualize your data. 
% Also, you can explore several helpful GUIs for data inspection
% and labeling.
%
% Created by Amirhossein Asadian
% Last Modified: Nov 29, 2024

% Here's a demo session (This is recording from SC)
clc; clear;

ses_name = 'Be231205';
cd(['C:\Users\CorneilLab\Desktop\AA\Belle\', ses_name])

region = 'SC';

load(['Z_', ses_name, '_', 'Analog.mat'])
load(['Z_', ses_name, '_', 'Digital.mat'])
load(['Z_', ses_name, '_', 'ML.mat'])
load(['Z_', ses_name, '_', region, '_Spikes.mat'])
load(['Z_', ses_name, '_', region, '_LFP.mat'])

folder_path = 'C:\Users\CorneilLab\Desktop\SortingTools';
addpath(genpath(folder_path));

folder_path = 'C:\Users\CorneilLab\Desktop\SharedUtils';
addpath(genpath(folder_path));

folder_path = 'C:\Users\CorneilLab\Documents\NIMH_MonkeyLogic_2.2';
addpath(genpath(folder_path)); 


%% Function name:    get_continuous_trial.m
% Purpose:          It gets a continous data and passes the aligned continous data
% Notes:            This is a general function for getting aligned
%                   continous data from your Z files. Read the help for more details.

tr_num = 1:100;
[cont_tr] = get_continuous_trial(Z_Analog, Z_Analog.data(1:2, :), tr_num, 'PDnumber', 1, 'fs', 1000, 'predur', 300, 'postdur', 1000);

plot(-300:1000, reshape(cont_tr(1, 1, :), 1, size(cont_tr, 3)))



%% Function name:    get_spikes_trial.m
% Purpose:          It gets spiking data and passes the aligned spikes and PSTH
% Notes:            This is a general function for getting aligned
%                   spiking and PSTH data from your Z files. Read the help for more details.

cond_want = 5;          %5: ETP
trs_right = find(Z_ML.condition == cond_want & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Right'));
ch = 3;
unitID = 1;
spk = Z_Spikes.data(ch).neural_timeStamps(1, find(Z_Spikes.data(ch).neural_timeStamps(2, :) == Z_Spikes.data(ch).Units(unitID)));
[SPIKES_RIGHT, SPDEN_RIGHT] = get_spikes_trial(Z_Analog, spk, trs_right, 'PDnumber', 1, 'predur', 300, 'postdur', 1000);

plot(-300:1000, mean(SPDEN_RIGHT))



%% Function name:    analyze_NHP_plotAll.m
% Purpose:          Plot spdn for all recorded neurons
% Notes:            Using this function you can plot PSTH for all neurons
%                   for left and rightward conditions. This function uses
%                   slider (something you might want to use in other
%                   places!). Alignment can be 'ReachGO' 'Target' 'Saccade'
%                   'Arm'. ch_want is a range of channels you want to plot.
%                   This function is designed for SLR project.

ch_want = 1:32;
cond_want = 5;
alignment = 'Target';
analyze_NHP_plotAll(Z_ML, Z_Analog, Z_Spikes, ch_want, cond_want, alignment)



%% Function name:    my_eyetraceViewer.m
% Purpose:          GUI for labeling eye trajectories
% Notes:            This GUI shows right and leftward trajectories with
%                   saccade RT. You can label the saccade movement. Check
%                   the function text right after you run it for key press
%                   behavior.

mystruct.Z_Name  = ['Z_', ses_name, '_ML'];

my_eyetraceViewer(Z_ML, mystruct);



%% Function name:    my_armtraceViewer.m
% Purpose:          GUI for labeling touchscreen movement trajectories
% Notes:            This GUI shows right and leftward trajectories with
%                   Arm RT. You can label the Arm movement. Check
%                   the function text right after you run it for key press
%                   behavior.

mystruct.Z_Name  = ['Z_', ses_name, '_ML'];

my_armtraceViewer(Z_ML, mystruct);



%% Function name:    get_SC_map.m
% Purpose:          Plot SC retinotropic map. Based on Ottes et al., 1986.

get_SC_map()



%% Function name:    my_psthViewer.m
% Purpose:          GUI for inspecting and labeling units recorded from the
%                   motor cortex
% Notes:            This GUI plots PSTH, spikes and waveforms for units
%                   recorded from the motor cortex. You can label units
%                   and inspect the ROCs. The receptive field, and
%                   classification tasks are also shown in this GUI.

cond_want = 5;          %5: ETP
trs_right = find(Z_ML.condition == cond_want & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Right'));
trs_left = find(Z_ML.condition == cond_want & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Left'));
% spikeTimes = Z_Spikes.data(ch).neural_timeStamps(1, find(Z_Spikes.data(ch).neural_timeStamps(2, :) == Z_Spikes.data(ch).Units(unitID)));

event_times_r = Z_Analog.PD_sec(trs_right);
event_times_l = Z_Analog.PD_sec(trs_left);
event_times_r_eye = event_times_r +  Z_ML.EyeRT(trs_right)/1000;
event_times_l_eye = event_times_l +  Z_ML.EyeRT(trs_left)/1000;
event_times_r_arm = event_times_r +  Z_ML.ArmRT(trs_right)/1000;
event_times_l_arm = event_times_l +  Z_ML.ArmRT(trs_left)/1000;

trs_eyeonly = find((Z_ML.condition == 1 | Z_ML.condition == 3) & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Right'));
trs_eyehand = find((Z_ML.condition == 2 | Z_ML.condition == 4) & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Right'));
event_times_eyeonly = Z_Analog.PD_sec(trs_eyeonly);
event_times_eyehand = Z_Analog.PD_sec(trs_eyehand);
event_times_eyeonly_eye = event_times_eyeonly +  Z_ML.EyeRT(trs_eyeonly)/1000;
event_times_eyehand_eye = event_times_eyehand +  Z_ML.EyeRT(trs_eyehand)/1000;

event_times = [event_times_r, event_times_l, event_times_r_eye, event_times_l_eye, ...
    event_times_r_arm, event_times_l_arm, event_times_eyeonly, event_times_eyehand, event_times_eyeonly_eye, event_times_eyehand_eye];

trialGroups = [zeros(size(event_times_r))+1, zeros(size(event_times_l))+2, ...
    zeros(size(event_times_r_eye))+3, zeros(size(event_times_l_eye))+4, ...
    zeros(size(event_times_r_arm))+5, zeros(size(event_times_l_arm))+6, ...
    zeros(size(event_times_eyeonly))+7, zeros(size(event_times_eyehand))+8, ...
    zeros(size(event_times_eyeonly_eye))+9, zeros(size(event_times_eyehand_eye))+10];

window = [-0.3 1]; % look at spike times from 0.3 sec before each event to 1 sec after
mystruct.Z_Name  = ['Z_', ses_name, '_', region, '_Spikes'];
my_psthViewer(Z_Spikes, Z_Analog, event_times, trialGroups, window, 1, 0, 1, mystruct);
% my_psthViewer(Z_Spikes, Z_Analog, eventTimes, trGroups, window, plotclassification, plotdelayed, plotmap, varargin)



%% Function name:    my_psthViewer_Frontal.m
% Purpose:          GUI for inspecting and labeling units recorded from the
%                   motor cortex
% Notes:            This GUI plots PSTH, spikes and waveforms for units
%                   recorded from the motor cortex. You can label units
%                   and inspect the ROCs. The receptive field, and
%                   classification tasks are not shown in this GUI compared
%                   to the SC one.

cond_want = 5;          %5: ETP
trs_right = find(Z_ML.condition == cond_want & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Right'));
trs_left = find(Z_ML.condition == cond_want & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Left'));
% spikeTimes = Z_Spikes.data(ch).neural_timeStamps(1, find(Z_Spikes.data(ch).neural_timeStamps(2, :) == Z_Spikes.data(ch).Units(unitID)));

event_times_r = Z_Analog.PD_sec(trs_right);
event_times_l = Z_Analog.PD_sec(trs_left);
event_times_r_eye = event_times_r +  Z_ML.EyeRT(trs_right)/1000;
event_times_l_eye = event_times_l +  Z_ML.EyeRT(trs_left)/1000;
event_times_r_arm = event_times_r +  Z_ML.ArmRT(trs_right)/1000;
event_times_l_arm = event_times_l +  Z_ML.ArmRT(trs_left)/1000;

trs_eyeonly = find((Z_ML.condition == 1 | Z_ML.condition == 3) & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Right'));
trs_eyehand = find((Z_ML.condition == 2 | Z_ML.condition == 4) & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Right'));
event_times_eyeonly = Z_Analog.PD_sec(trs_eyeonly);
event_times_eyehand = Z_Analog.PD_sec(trs_eyehand);
event_times_eyeonly_eye = event_times_eyeonly +  Z_ML.EyeRT(trs_eyeonly)/1000;
event_times_eyehand_eye = event_times_eyehand +  Z_ML.EyeRT(trs_eyehand)/1000;

event_times = [event_times_r, event_times_l, event_times_r_eye, event_times_l_eye, ...
    event_times_r_arm, event_times_l_arm, event_times_eyeonly, event_times_eyehand, event_times_eyeonly_eye, event_times_eyehand_eye];

trialGroups = [zeros(size(event_times_r))+1, zeros(size(event_times_l))+2, ...
    zeros(size(event_times_r_eye))+3, zeros(size(event_times_l_eye))+4, ...
    zeros(size(event_times_r_arm))+5, zeros(size(event_times_l_arm))+6, ...
    zeros(size(event_times_eyeonly))+7, zeros(size(event_times_eyehand))+8, ...
    zeros(size(event_times_eyeonly_eye))+9, zeros(size(event_times_eyehand_eye))+10];

window = [-0.3 1]; % look at spike times from 0.3 sec before each event to 1 sec after
mystruct.Z_Name  = ['Z_', ses_name, '_', region, '_Spikes'];
% my_psthViewer(Z_Spikes, Z_Analog, event_times, trialGroups, window, 0, 0, 0, mystruct);
my_psthViewer_Frontal(Z_Spikes, Z_Analog, event_times, trialGroups, window, 0, 0, 0, mystruct)
% my_psthViewer(Z_Spikes, Z_Analog, eventTimes, trGroups, window, plotclassification, plotdelayed, plotmap, varargin)


%% Function name:    my_csdViewer.m
% Purpose:          plot CSD for LFP
% Notes:            This function shows CSD for a matrix of LFPs over
%                   channels. Also, it gives you the csd profiles of each
%                   trials and averages. It is based on the CSDplotter
%                   toolbox.

[cont_tr] = get_continuous_trial(Z_Analog, Z_LFP.data(:, :), find(Z_ML.TrialError == 1 & Z_ML.condition == 5));

[b, a] = butter(3, [1, 250]/(Z_Analog.info.SampleRate/2), 'bandpass'); % lowpass filter
    
[b2, a2] = butter(3, [0.7, 30]/(Z_Analog.info.SampleRate/2), 'bandpass'); % lowpass filter

% Bandstop (notch) filter for 60 Hz line noise removal
d = designfilt('bandstopiir','FilterOrder',6, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',Z_Analog.info.SampleRate);

for i = 1:size(cont_tr, 1)
    cont_tr_ch = reshape(cont_tr(i, :, :), size(cont_tr, 2), size(cont_tr, 3));
    cont_tr_ch = filtfilt(b, a, cont_tr_ch');
    cont_tr_ch = filtfilt(d, cont_tr_ch);
    cont_tr_ch = filtfilt(b2, a2, cont_tr_ch);
    
    cont_tr(i, :, :) = cont_tr_ch';
end

my_csdViewer(cont_tr, 50, 150)


%% Function name:    disp_LFP_inline.m
% Purpose:          plot LFP over channels
% Notes:            Simple plot function, plotting a continous data over
%                   channels

cond_want = 5;          %5: ETP
trs_right = find(Z_ML.condition == cond_want & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Right') & [Z_ML.EyeRT(:, 1); zeros(length(Z_ML.TrialError) - length(Z_ML.EyeRT(:, 1)), 1)]' >= 60);
trs_left = find(Z_ML.condition == cond_want & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Left') & [Z_ML.EyeRT(:, 1); zeros(length(Z_ML.TrialError) - length(Z_ML.EyeRT(:, 1)), 1)]' >= 60);

[cont_tr_r] = get_continuous_trial(Z_Analog, Z_LFP.data(:, :), trs_right, 'predur', 500, 'postdur', 1000);
[cont_tr_l] = get_continuous_trial(Z_Analog, Z_LFP.data(:, :), trs_left, 'predur', 500, 'postdur', 1000);

[b, a] = butter(3, [1, 250]/(Z_Analog.info.SampleRate/2), 'bandpass'); % lowpass filter
    
[b2, a2] = butter(3, [0.7, 30]/(Z_Analog.info.SampleRate/2), 'bandpass'); % lowpass filter

% Bandstop (notch) filter for 60 Hz line noise removal
d = designfilt('bandstopiir','FilterOrder',6, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',Z_Analog.info.SampleRate);

for i = 1:size(cont_tr_r, 1)
    cont_tr_ch = reshape(cont_tr_r(i, :, :), size(cont_tr_r, 2), size(cont_tr_r, 3));
    cont_tr_ch = filtfilt(b, a, cont_tr_ch');
    cont_tr_ch = filtfilt(d, cont_tr_ch);
    cont_tr_ch = filtfilt(b2, a2, cont_tr_ch);
    cont_tr_r(i, :, :) = cont_tr_ch';
end

for i = 1:size(cont_tr_l, 1)
    cont_tr_ch = reshape(cont_tr_l(i, :, :), size(cont_tr_l, 2), size(cont_tr_l, 3));
    cont_tr_ch = filtfilt(b, a, cont_tr_ch');
    cont_tr_ch = filtfilt(d, cont_tr_ch);
    cont_tr_ch = filtfilt(b2, a2, cont_tr_ch);
    cont_tr_l(i, :, :) = cont_tr_ch';
end

lfp_data_r = reshape(mean(cont_tr_r, 2), size(cont_tr_r, 1), []);
lfp_data_l = reshape(mean(cont_tr_l, 2), size(cont_tr_l, 1), []);

figure
disp_LFP_inline(gca, lfp_data_r, 500, 1000, 100, 1, [], '', 'b');
disp_LFP_inline(gca, lfp_data_l, 500, 1000, 100, 1, [], '', 'r');
