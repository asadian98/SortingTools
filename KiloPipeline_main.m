%% Data extraction and spike sorting pipeline, CorneilLab, Western University
% Created by Amirhossein Asadian
% Last Modified: Dec 12, 2024

% Notes:
% This is a pipeline script for extracting data from MonkeyLogic and Ripple,
% as well as sorting neural data using Kilosort2. You can extract
% various files based on your selection in the following lines.

% **** Raw neural data files should be in this directory:
% Desktop --> {Researcher} --> Ripple_Raw_Data --> {animalName} --> {sesName}{tail}.nev/.ns2/.ns5
% **** MonkeyLogic files should be in this directory:
% Desktop --> {Researcher} --> ML_Behavioural_Files --> {animalName} --> {sesName}{tail}.bhv2

% **** Output files will be in this directory:
% Desktop --> {Researcher} --> {animalName} --> {sesName} --> (outputfiles)

% outputfiles:
%             Z_{session name}_Analog.mat
%                                           Z_Analog:
%                                                       data (raw data for eye_x, eye_y, eye_pupil, photodiode) --> dimension is 4*numsamples
%                                                       info
%                                                       PD_sec (photodiode onset times relative to the beginning of the session) --> numtrials*numPDs
%
%             Z_{session name}_Digital.mat
%                                           Z_Digital:
%                                                       data (statecodes; statecodetimes) --> dimension is 2*numstatecodes
%                                                       info
%
%             Z_{session name}_ML.mat
%                                           Z_ML (MonkeyLogic data including calibrated eye, touch, ... )
%
%             Z_{session name}_{region name}_LFP.mat
%                                           Z_LFP
%
%             Z_{session name}_{region name}_Spikes.mat
%                                           Z_Spikes
%                                                       data: spikes and times of spikes for each channel (the contact on which units sit comes from Kilosort)
%                                                       info
%                                                       su: information for each unit (note that clusterId starts from zero, the peakCh might be different than the channel number in data)
%
%             Z_{session name}_{region name}_TR.mat
%                                           Z_TR
%                                                       data: DLC tracking data for each frame
%                                                       frames_timestamps: frame timestamps
%                                                       ledvalue: pixel brightness related to photodiode on each frame
% You can also have Z_Wave/Z_Raw (very large files!)

% TODO:
%       - Bombcell

clc; clear;

folder_path = 'C:\Users\CorneilLab\Desktop\SortingTools';
addpath(genpath(folder_path));

folder_path = 'C:\Users\CorneilLab\Desktop\SharedUtils';
addpath(genpath(folder_path));

folder_path = 'C:\Users\CorneilLab\Documents\NIMH_MonkeyLogic_2.2';
addpath(genpath(folder_path)); 

sesName =                   'Be231201';
region  =                   'SC';                      % PMd, M1 or SC
tail =                      '_001';                     % session names should be like: Be231124_002;
% If you are partitioning ML or ripple file and want to concatanate them
% post recording, you need to change the ML_list or tail list

location =                  [4, 4];                     % posterior/Lateral
ch_offset =                 0;                          % If you recorded multiple electrode, the second electrode starts from 32 (A1-A2-...)

animalName =                'Belle';
Researcher =                'AA';

S_probe16 =                 0;                          % If 16 channels S-probe

getRPsum_flag =             0;                          % gives you ripple summary
getEvents_flag =            0;                          % get events
saveRaw_flag =              0;                          % Save raw neural ripple file
getDat_flag =               0;                          % save dat file
doFilters_flag =            0;                          % Apply filters
runKilo_flag =              1;                          % Run Kilosort
runQualityControl_flag =    0;                          % Run Bombcell for Quality Control
runPhy_flag =               1;                          % Run Phy after sorting
saveSpikes_flag =           0;                          % Save Z_Spikes (It also extract waveforms) extract waveform, quality metrics and ...; waveform extraction takes some time and will create the Z_WFf file
saveWave_flag =             0;                          % Save Z_Wave,
do_posthoc_flag =           0;                          % Just nice plots
get_DLC_flag =              0;                          % If you had behavioral tracking
KL_version =                4;                          % Kilosort version; 2 or 4

fs =                        30000;

% Specifiy and Extract Variables from Electrophysiology file
% Which file do you want to open and which electrodes are good?
Physio_Info.Statecodes =                1;
Physio_Info.Acute_EMG_DATA =            0;
Physio_Info.Good_AcuteEMG_Electrodes =  [0 0 0];
Physio_Info.GoodDay =                   1;
Physio_Info.RippleDataFile =            1;
Physio_Info.Pupil_Analog =              1;
Physio_Info.PlexonDataFile =            0;
Physio_Info.IntanDataFile =             0;
Physio_Info.NHPName =                   animalName;
Physio_Info.Analog_DATA =               1;
Physio_Info.AnalogSampleRate =          1000;
Physio_Info.Eye_Analog =                1;
Physio_Info.PD_Analog =                 1;
Physio_Info.R32_EMG_DATA =              0;
Physio_Info.Neural_DATA =               1;
Physio_Info.PostHocSpikeSort =          [];
Physio_Info.FHC_DATA =                  0;
Physio_Info.Plexon_Data =               1;
Physio_Info.Neural_Array1_Data =        0;
Physio_Info.Neural_Array2_Data =        0;
Physio_Info.PlotData =                  0;
Physio_Info.Sorter =                    'Kilosort';
Physio_Info.predur_min =                200;
Physio_Info.predur_max =                500;
Physio_Info.postdur =                   500;

if(S_probe16) nCh = 16; Physio_Info.numCh = 16;
else          nCh = 32; Physio_Info.numCh = 32;
end

ML_tail = tail;
Physio_Info.Concatenated = 0;
Physio_Info.DataFileName = sesName; % Physiology data file name

fileName = [sesName, tail];

Folders.O_NHP_Folder = ['C:\Users\CorneilLab\Desktop\', Researcher, '\'];
Folders.NHP_Folder = animalName;
Folders.Ripple_Raw_Data_Folder = [Folders.O_NHP_Folder,'Ripple_Raw_Data\'];
Folders.sortingFolder = [Folders.O_NHP_Folder, Folders.NHP_Folder, '\', sesName, '\', region, '_KS', num2str(KL_version)];
Folders.ML_Behav_Folder = [Folders.O_NHP_Folder,'\ML_Behavioural_Files\', animalName];
Folders.KiloFolder = [Folders.O_NHP_Folder, Folders.NHP_Folder, '\'];
Folders.save_dir = [Folders.O_NHP_Folder, Folders.NHP_Folder, '\', sesName, '\'];
Folders.TrackingUtils = ['C:\Users\CorneilLab\Desktop\SortingTools\Utils\Tracking\'];
Folders.Tracking = [Folders.O_NHP_Folder, 'Tracking\', animalName, '\'];
Folders.Utils = ['C:\Users\CorneilLab\Desktop\SortingTools\Utils\'];

cd(Folders.KiloFolder)
if ~exist([Folders.KiloFolder, sesName], 'dir')
    mkdir(sesName);
    disp(['Folder "' sesName '" was created.']);
else
    disp(['Folder "' sesName '" already exists.']);
end

cd([Folders.O_NHP_Folder, Folders.NHP_Folder, '\', sesName])
if ~exist([region, '_KS', num2str(KL_version)], 'dir')
    mkdir([region, '_KS', num2str(KL_version)]);
    disp(['Folder "' [region, '_KS', num2str(KL_version)] '" was created.']);
else
    disp(['Folder "' [region, '_KS', num2str(KL_version)] '" already exists.']);
end

% meta info
meta.nCh = nCh;
meta.S_probe16 = S_probe16;
meta.fs = fs;
meta.sesName = sesName;
meta.fileName = fileName;
meta.location = location;
meta.region = region;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      Ripple Summary     %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(getRPsum_flag)
    summary_ripple_
    recording(animalName, Researcher, sesName, tail)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Extract MonkeyLogic Data   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(getEvents_flag)

    % Here I'm creating a map (dict in python) defining different
    % tasks and conditions -- change this based on your MonkeyLogic script

    %             % For Grover
    %             keys = [502, 504, 121, 4];  % keys are monkeylogic codes
    %             values = struct('condition', {3, 4, 5, 1}, 'description', {'EYEHAND_EyeT1', 'EYEHAND_EyeT1_HandT2', 'Emerging Target: Single Target Red', 'EYE_EyeT1'});

    % For Belle
    keys = [500, 502, 121, 122, 620, 700, 120, 600];  % keys are monkeylogic codes
    values = struct('condition', {1, 2, 5, 6, 7, 8, 9, 10}, 'description', ...
        {'EYE_EyeT1', 'HAND_HandT2', 'Emerging Target: Single Target Red', 'Emerging Target: EyeOnly', 'MAPPING', 'Delayed Task', 'Emerging Target: static', 'center-out reaching'});

    % Initialize the map
    conditionCodeMap = containers.Map('KeyType', 'double', 'ValueType', 'any');

    for i = 1:length(keys)
        conditionCodeMap(keys(i)) = values(i);  % Assign struct as value for each key
    end

    % For Belle (For Grover GAP start and end are 40 and 30 respectively
    keys = 1:13;  % keys are monkeylogic codes
    values = struct('code', {9, 18, 30, 30, 60, 200, 200, 120, 110, 100, 124, 80, 35}, 'description', ...
        {'Trial Start', 'Trial End', 'Go Cue: EYEHAND_EyeT1', 'Go Cue: EYEHAND_EyeT1_HandT2', 'Emerging Target Task: Single Target', ...
        'Emerging Target Task: Double Target Red', 'Emerging Target Task: Double Target Green', 'Reach Go Cue: Delayed Task', 'Saccade Go Cue: Delayed Task', ...
        'Visual Go Cue: Delayed Task', 'Photodiode off', 'GAP start', 'GAP end'});

    % Initialize the map
    strobecodesMap = containers.Map('KeyType', 'double', 'ValueType', 'any');

    for i = 1:length(keys)
        strobecodesMap(keys(i)) = values(i);  % Assign struct as value for each key
    end

    Z = Concat_ML_Files(sesName, Physio_Info, {ML_tail}, {tail}, animalName, Folders, conditionCodeMap, strobecodesMap);

    % get Event times (you need this for psth viewer in PHY)
    getEvents_concat_segment(sesName, Z, animalName, {tail}, Physio_Info, Folders)

end

% --------------------------- Dat File --------------------------------
if(getDat_flag)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%     Extract Z Files     %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Z = [];
    Z_Raw = Concat_Ripple_Extract_Signals(Z, Physio_Info, 32, sesName, {tail}, ch_offset, region, Folders, saveRaw_flag);

    RawNeural_Data_mV = Z_Raw.data;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%    Filtering    %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(doFilters_flag)
        RawNeural_Data_mV = doFilters(RawNeural_Data_mV, fs);
    end

    % ---------------------- Save Dat File ----------------------------

    cd(Folders.KiloFolder)
    mkdir('temp');
    datFile = [Folders.KiloFolder, sesName, '\', fileName, '.dat'];

    % write ephys data to dat file:
    fidout = fopen(datFile, 'w'); % opening file for appending
    fwrite(fidout, RawNeural_Data_mV, 'int16');

    fclose(fidout);
    disp('Conversion complete!')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Run Kilosort (v2 and v4)  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(runKilo_flag)
    if(KL_version == 2)
        runKilo(Folders, meta, -6); % -6 is the default threshold
    elseif(KL_version == 4) % You need to find 'good' parameters for your recording at first
        cd(Folders.Utils)
        anaconda_path = 'C:\ProgramData\anaconda3';
        setenv('PATH', [anaconda_path ';' anaconda_path '\Scripts;' anaconda_path '\Library\bin;' getenv('PATH')]);
        if(~meta.S_probe16)
            command = strcat("activate kilosort && python run_KS4.py ", Folders.KiloFolder, sesName, "\ ", Folders.sortingFolder, " ", Folders.Utils, "chanMap_32ch_150um.mat ", "32 ", meta.fileName,".dat");
        else
            command = strcat("activate kilosort && python run_KS4.py ", Folders.KiloFolder, sesName, "\ ", Folders.sortingFolder, " ", Folders.Utils, "chanMap_16ch_300um.mat ", "16 ", meta.fileName,".dat");
        end
        [status, cmdout] = system(command, "-echo")
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Quality Control - Bombcell  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% right now, I'm not using any additional quality control. getSp function
% in saveSpikes script calculates violations. There are more quality control
% in Bombcell toolbox. There are also some useful waveform extraction
% functions in that toolbox as well.

if(runQualityControl_flag)
    runQualityControl(Folders, meta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%            PHY            %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Additional features we've added to our PHY:
% Rasters and PSTH plotting for different conditions and events
% Outlier removal using the Mahalanobis distance. Standard threshold is 16 standard deviations (adjustable threshold).
% Reclustering with KlustaKwik 2.0 and kmeans

if(runPhy_flag)
    anaconda_path = 'C:\ProgramData\anaconda3';
    setenv('PATH', [anaconda_path ';' anaconda_path '\Scripts;' anaconda_path '\Library\bin;' getenv('PATH')]);
    command = strcat("activate phy2 && phy template-gui  ", Folders.sortingFolder, "\params.py");
    [status, cmdout] = system(command, "-echo")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get and Save Spikes from PHY  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(saveSpikes_flag)
    saveSpikes(Folders, meta, Physio_Info, saveWave_flag);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Post-hoc Visualization    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(do_posthoc_flag)
    do_posthoc(Folders, meta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Behavioural Tracking     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Takes about 45 mins to run DLC
if(get_DLC_flag)
    getDLC(Folders, meta);
end
