%% Ripple Physiology Data Extraction
% Created by Aaron L. Cecala (ALC), modified by Amirhossein Asadian 
% Last Modified: August 22, 2024

% This function does not extract spikes for you, spikes come from Kilosort.

% Notes:
% This is a function template for accessesing different data from a Trellis
% data set and uses a specific version of the NeuroShare toolbox that comes with the
% Trellis software package. Will not work otherwise. This is called
% "Corneil_Neuroshare" in ALC's original file structure.

% If you are using a data set that utilizes pausing, the time between data
% will be padded with ZEROS. However, this is typically not how the Corneil lab
% collects data.

% Output from this function will be both a "raw" datastream from the
% representing all the data from the '*.NS#' and '.NEV' datafiles from
% a particular session as well as spliced data that represents data
% from each trial as based on strobe codes. Further processing of data
% (e.g., alignment of data to a particular event) must occur in another
% function.

% Inputs:
% O_NHP_Folder - path to organize folder

% Physio_Info - a structure that contains:
% RippleDataFile = 1;  1 --> this is a Ripple Datafile
% NHPName = 'Grover';  Name of non-human primate subject
% DataFileName = '06242022E'; '*.NS#' datafile string
% Statecodes = 1; 1 --> you want statecodes extracted
% Analog_DATA = 1; 1 --> you want analog data extracted
% AnalogSampleRate = 1000; 1000 = 1K; 30000 = 30K sample rate
% Eye_Analog = 1; 1 --> you want uncalibrated eye extracted
% PD_Analog = 1; 1 --> you want photodiode data extracted
% Acute_EMG_DATA = 0; 1 --> you want acute EMG extracted
% Good_AcuteEMG_Electrodes = [0 0 0]; This is an array of zeros
% and ones where one implies the electrode was a good electrode
% for data processing and zero is a bad electrode. Length
% of array must equal number of "muscles' sampled.
% R32_EMG_DATA = 0; 1 --> you want potted link R32 data extracted
% Neural_DATA = 1; 1 --> you want potted Neural data extracted
% FHC_DATA = 0; 1 --> you want single contact electrode data
% processed (currently not written in this code, but would be
% first channel from omnetics connector.)
% Plexon_Data = 1; 1 --> you want to process S/V-Probe ephys
% Neural_Array1_Data = 0; 1 --> you want to process chronic
% neural ephys from array 1. (currently not written in this
% code)
% Neural_Array2_Data = 0; 1 --> you want to process chronic
% neural ephys from array 2. (currently not written in this
% code)
% PlotData = 0; 1 --> you want to plot anything as a test of
% this script.

% OUTPUTS:
% PHYSIOLOGY_RAW - contains full file stream of selected
% trial information and physiology data.

function Z_Raw = Concat_Ripple_Extract_Signals(Z, Physio_Info, numberOfChan, ses_name_main, tail_list, ch_offset, region, Folders, saveRaw)

% Concat has not been modified for EMGs
disp('EXTRACTING RIPPLE SIGNALS')
time_offset = 0;

if(strcmp(Physio_Info.NHPName, 'Grover') && ~strcmp(ses_name_main(1), 'G')) % Grover's data channel numbers starts from 257
    ch_offset = 256; % 256 for Aaron's recording
end

for tail_idx = 1:length(tail_list)

    ses_name = [ses_name_main, tail_list{tail_idx}];


    % What RIPPLE data do you want to look at?
    DataFileName = ses_name;
    Statecodes = Physio_Info.Statecodes;
    AnalogSampleRate = Physio_Info.AnalogSampleRate;

    Analog_Data = Physio_Info.Analog_DATA;
    Eye_Analog = Physio_Info.Eye_Analog;
    PD_Analog = Physio_Info.PD_Analog;

    Acute_EMG_DATA = Physio_Info.Acute_EMG_DATA;
    Good_AcuteEMG_Electrodes = Physio_Info.Good_AcuteEMG_Electrodes;
    R32_EMG_DATA = Physio_Info.R32_EMG_DATA;

    Neural_DATA = Physio_Info.Neural_DATA;
    FHC_DATA = Physio_Info.FHC_DATA;
    Plexon_Data = Physio_Info.Plexon_Data;
    Neural_Array1_Data = Physio_Info.Neural_Array1_Data;
    Neural_Array2_Data = Physio_Info.Neural_Array2_Data;

    % Conversion Factors (this is to make sure that Neuroshare output is in
    % microvolts for PottedLink 32)
    units_multiplier = 5; % Checked with Blackrock toolbox used by Kevin Cross (Queens).

    %% Data Path and FileName Information for RIPPLE FILE
    % Data File path
    DataPath = 'Ripple_Raw_Data\';
    DataFolder = Physio_Info.NHPName;

    %% Get Statecodes from Digital Datafile
    if Statecodes == 1
        disp('Getting Digital Statecodes for Behavioural Analysis')
        % Complete Data File Path: Digital Data
        completeFilePath_Digital = [Folders.O_NHP_Folder,DataPath,DataFolder,'\',DataFileName,'.nev']
        % Open the file and extract some basic information
        disp('opening Digital datafile to acquire statecodes')
        [ns_status, hFile] = ns_OpenFile(completeFilePath_Digital,'single');  % opens only .nev file
        if strcmp(ns_status,'ns_OK')
            disp('Digital datafile opened properly')
        else
            disp('Did not open Digital datafile properly')
        end
        % ElectrodeID of ZERO is the Digital Input Parallel port input
        EntityIndices = find([hFile.Entity(:).ElectrodeID] == 0);
        entityID = EntityIndices(1);
        [ns_RESULT, entityInfo] = ns_GetEntityInfo(hFile, entityID(end));

        % Get events and time stamps

        % TimeStamp           Variable that receives the timestamp of the Event
        %                     data item.
        %
        % Data                Variable that receives the data for the Event entry.
        %                     The format of data is specified by the member
        %                     EventType in ns_EVENTINFO.
        %
        % DataSize            Variable that receives the actual number of bytes of
        %                     data retrieved in the data buffer.

        numCount = entityInfo.ItemCount; % how many statecode events were there?
        % create dummy arrays to fill with statecode information
        statecodes = NaN(1, numCount); statecode_timeStamps = NaN(1, numCount); dataSize = NaN(1, numCount);
        for i = 1:numCount
            [~, statecode_timeStamps(i), statecodes(i), dataSize(i)] = ns_GetEventData(hFile, entityID, i);
        end

        % Find Unique Statecodes (sanity check)
        unique_statecodes = unique(statecodes);
        if ~isempty(unique_statecodes)
            disp('Statecodes were extracted properly')
            disp('Statecodes:')
            for statecode_num = 1:length(unique_statecodes)
                disp(num2str(unique_statecodes(statecode_num)))
            end
        else
            disp('PROBLEM: There were no unique statecodes!')
        end

        if(statecodes(1) == statecodes(2))
            statecodes = statecodes(1:2:end);
            statecode_timeStamps = statecode_timeStamps(1:2:end);
        end

        num9 = find(statecodes == 9);
        num18 = find(statecodes == 18);

        % You should comment this part if you have segmented recording
%         if(num9(1) ~= 1)
%             statecodes(1:num9(1)-1) = [];
%             statecode_timeStamps(1:num9(1)-1) = [];
%         end
% 
%         if(num18(end) ~= length(statecodes))
%             statecodes((num18(end) + 1) : end) = [];
%             statecode_timeStamps((num18(end) + 1) : end) = [];
%         end
% 
%         % For files that have less trials in ML file
%         st_idx = find(statecodes == 9);
%         st_idx2 = find(statecodes == 18);
%         if(length(st_idx) > Z.tr_sep(tail_idx))
%             statecodes(st_idx2(Z.tr_sep(tail_idx))+1:end) = [];
%             statecode_timeStamps(st_idx2(Z.tr_sep(tail_idx))+1:end) = [];
%         end
%         assert(length(find(statecodes == 9)) == Z.tr_sep(tail_idx), 'Error: check number of trials, somethings wrong');


        if(tail_idx == 1)
            DigitalInfo.Statecodes = [];
            DigitalInfo.Statecodes = statecodes;
            DigitalInfo.Statecode_timeStamps = statecode_timeStamps;
            %             DigitalInfo.DataSize = dataSize;
        else
            DigitalInfo.Statecodes = [DigitalInfo.Statecodes, statecodes];
            DigitalInfo.Statecode_timeStamps = [DigitalInfo.Statecode_timeStamps, statecode_timeStamps + time_offset];
            %             DigitalInfo.DataSize = [DigitalInfo.DataSize, dataSize];
        end
    else
        disp('Did NOT get Digital Statecodes for Behavioural Analysis')
    end

    %% Extract Analog Inputs (EyeX, EyeY, Pupil, Photodiode(PD), Acute EMGs)
    if Analog_Data == 1
        disp('Extracting Analog Channels for Behavioural Analysis')
        if AnalogSampleRate == 30000
            % Complete Data File Path
            completeFilePath_Analog = strcat(Folders.O_NHP_Folder,DataPath,DataFolder,'/',DataFileName,'.ns5');
            dataType_ANALOG = 'Analog 30k';      % could also be 'Analog 1k' depending on Trelis data collect settings
            Info.SampleRate = 30000;
        elseif AnalogSampleRate == 1000
            completeFilePath_Analog = strcat(Folders.O_NHP_Folder,DataPath,DataFolder,'/',DataFileName,'.ns2');
            dataType_ANALOG = 'Analog 1k';      % could also be 'Analog 30k' depending on Trelis data collect settings
            Info.SampleRate = 1000;
        end

        disp(dataType_ANALOG)

        Analog_dataChannels = [10241, 10242, 10243, 10244];
        % Analog I/O starts at 10241;
        % 10241 = horizontal eye raw, uncalibrated
        % 10242 = vertical eye raw, uncalibrated
        % 10243 = pupil raw, uncalibrated
        % 10244 = photodiode
        % 10245 = Acute EMG 1
        % 10246 = Acute EMG 2
        % 10247 = Acute EMG 3

        numAnalog_Channels = length(Analog_dataChannels);

        %% Open the analog file and extract some basic information
        disp('opening analog datafile')
        [ns_status, hFile] = ns_OpenFile(completeFilePath_Analog,'single');
        if strcmp(ns_status,'ns_OK')
            disp('Analog datafile opened properly')
        else
            disp('Did not open Analog datafile properly')
        end

        EntityIndices = [];
        fileType = {};
        entityID = [];
        fileTypeNum = [];

        %% Determine correct entityID for desired datastream
        for analog_channel_num = 1:numAnalog_Channels
            EntityIndices(analog_channel_num) = find([hFile.Entity(:).ElectrodeID] == Analog_dataChannels(analog_channel_num));
            fileTypeNum(analog_channel_num) = hFile.Entity(EntityIndices(analog_channel_num)).FileType;
            fileType{analog_channel_num} = hFile.FileInfo(1).Type;
            entityID(analog_channel_num) = EntityIndices(analog_channel_num);
            % Extract channel info
            [ns_RESULT, entityInfo(analog_channel_num)] = ns_GetEntityInfo(hFile, entityID(end));
        end

        %% Extract all 30k or 1k analog data channels (in basic Ripple setup, there are only 4)
        disp('Extracting Analog Data')
        if(tail_idx ~= 1)
            tmp = Analog_Data_Time_s;
            tmp2 = Analog_Data_mV;
        end
        Analog_Data_Time_s = [];
        Analog_Data_mV = [];
        Analog_Data_Time_s2 = [];
        Analog_Data_mV2 = [];
        %         for channel_num_AnalogData = 1:length(EntityIndices)
        %             [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, entityID(channel_num_AnalogData));     % analog info contains things like range and sampling rate
        %
        %             TimeStamps = hFile.FileInfo(hFile.Entity(entityID(channel_num_AnalogData)).FileType).TimeStamps;
        %             numSamples = sum(TimeStamps(:,end));
        %             analogInputData = zeros(1,numSamples);
        %             startIndex = 1;
        %             indexCount = TimeStamps(2,1);
        %             for i = 1:size(TimeStamps,2)
        %                 [~, ~, tempData] = ns_GetAnalogData(hFile, entityID(channel_num_AnalogData), startIndex, indexCount);
        %                 dataRange = TimeStamps(1,i) + (1:TimeStamps(2,i));
        %                 analogInputData(dataRange) = tempData';
        %                 clear tempData
        %                 if i ~= size(TimeStamps,2)
        %                     startIndex = startIndex + TimeStamps(2,i);
        %                     indexCount = TimeStamps(2,i+1);
        %                 end
        %             end
        %             analogInputDataTime_s = (0:numSamples-1)' ./ analogInfo.SampleRate;
        %             analogInputData_mV = analogInputData; clear analogInputData;
        %
        %             if(tail_idx == 1)
        %                 Analog_Data_Time_s(channel_num_AnalogData,:) = analogInputDataTime_s';
        %                 Analog_Data_mV(channel_num_AnalogData,:) = analogInputData_mV;
        %             else
        %                 Analog_Data_Time_s(channel_num_AnalogData,:) = [tmp(channel_num_AnalogData, :), analogInputDataTime_s' + time_offset];
        %                 Analog_Data_mV(channel_num_AnalogData,:) = [tmp2(channel_num_AnalogData, :), analogInputData_mV];
        %             end
        %
        %         end

        analoginfo.scale = [hFile.Entity(EntityIndices(:)).Scale];
        for entIdx = 1:length(EntityIndices)
            analoginfo.units{entIdx} = [hFile.Entity(EntityIndices(entIdx)).Units];
        end

        if(tail_idx == 1)
            nTimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices(1)).FileType).TimeStamps;
            [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices(1));
            [ns_RESULT, Analog_Data_mV_tmp] = ns_GetAnalogDataBlock(hFile, EntityIndices, 1, nTimeStamps(2,1), 0);
            Analog_Data_Time_s(1:length(EntityIndices),:) = repmat((0:nTimeStamps(2,1)-1) ./ analogInfo.SampleRate, [length(EntityIndices), 1]);
            Analog_Data_mV(1:length(EntityIndices),:) = Analog_Data_mV_tmp';
        else
            nTimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices(1)).FileType).TimeStamps;
            [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices(1));
            [ns_RESULT, Analog_Data_mV_tmp] = ns_GetAnalogDataBlock(hFile, EntityIndices, 1, nTimeStamps(2,1), 0);
            Analog_Data_Time_s2(1:length(EntityIndices),:) = repmat((0:nTimeStamps(2,1)-1) ./ analogInfo.SampleRate, [length(EntityIndices), 1]);
            Analog_Data_mV2(1:length(EntityIndices),:) = Analog_Data_mV_tmp';
            Analog_Data_Time_s(1:length(EntityIndices),:) = [tmp(:, :), Analog_Data_Time_s2 + time_offset];
            Analog_Data_mV(1:length(EntityIndices),:) = [tmp2(:, :), Analog_Data_mV2];
        end

%         Eye_ArrayData_Stream = {}; % Is this correct? no tail?
%         if Eye_Analog == 1
%             Eye_ArrayData_Stream{1,1}(1,:) = Analog_Data_mV(1,:); % Horizontal Eye data in mV
%             Eye_ArrayData_Stream{1,1}(2,:) = Analog_Data_Time_s(1,:); % timestamps (depends on 1K vs. 30K)
%             Info.EyeChannel{1,1} = 'EyeX';
%             Eye_ArrayData_Stream{2,1}(1,:) = Analog_Data_mV(2,:); % Vertical Eye data in mV
%             Eye_ArrayData_Stream{2,1}(2,:) = Analog_Data_Time_s(2,:); % timestamps (depends on 1K vs. 30K)
%             Info.EyeChannel{2,1} = 'EyeY';
%             Eye_ArrayData_Stream{3,1}(1,:) = Analog_Data_mV(3,:); % Pupil data in mV
%             Eye_ArrayData_Stream{3,1}(2,:) = Analog_Data_Time_s(3,:); % timestamps (depends on 1K vs. 30K)
%             Info.EyeChannel{3,1} = 'Pupil';
%         else
%             Info.EyeChannel{1,1} = 'EyeX Empty';
%             Info.EyeChannel{2,1} = 'EyeY Empty';
%             Info.EyeChannel{3,1} = 'PupilEmpty';
%         end
% 
%         PD_ArrayData_Stream = {};
%         if PD_Analog == 1
%             PD_ArrayData_Stream{1,1}(1,:) = Analog_Data_mV(4,:); % Photodiode data in mV
%             PD_ArrayData_Stream{1,1}(2,:) = Analog_Data_Time_s(4,:); % timestamps (depends on 1K vs. 30K)
%             Info.PDChannel{1,1} = 'PhotoDiode Aligned to First Peripheral Target Onset';
%         end

        % _________________________________________________________________
        % ------------------------------ EMG ------------------------------
        % _________________________________________________________________
         
        emginfo = [];
        if(Acute_EMG_DATA)
            disp('Extracting 3k Analog Channels')
            completeFilePath_Analog = strcat(Folders.O_NHP_Folder,DataPath,DataFolder,'/',DataFileName,'.ns5');
            dataType_ANALOG = 'Analog 30k';      % could also be 'Analog 1k' depending on Trelis data collect settings
            Info2.SampleRate = 30000;

            disp(dataType_ANALOG)

            Analog_dataChannels = [10245, 10246, 10247];

            % Analog I/O starts at 10241;
            % 10241 = horizontal eye raw, uncalibrated
            % 10242 = vertical eye raw, uncalibrated
            % 10243 = pupil raw, uncalibrated
            % 10244 = photodiode
            % 10245 = Acute EMG 1
            % 10246 = Acute EMG 2
            % 10247 = Acute EMG 3

            numAnalog_Channels = length(Analog_dataChannels);

            %% Open the analog file and extract some basic information
            disp('opening analog datafile')
            [ns_status, hFile] = ns_OpenFile(completeFilePath_Analog,'single');
            if strcmp(ns_status,'ns_OK')
                disp('Analog datafile opened properly')
            else
                disp('Did not open Analog datafile properly')
            end

            EntityIndices = [];
            fileType = {};
            entityID = [];
            fileTypeNum = [];

            %% Determine correct entityID for desired datastream
            
            EntityIdx = 0;
            EntityIdx2 = 0;

            Entitynull = zeros(1, numAnalog_Channels);
            for analog_channel_num = 1:numAnalog_Channels
                EntityIdx2 = EntityIdx2 + 1;

                if(~isempty(find([hFile.Entity(:).ElectrodeID] == Analog_dataChannels(analog_channel_num))))
                    EntityIdx = EntityIdx + 1;
                    EntityIndices(EntityIdx) = find([hFile.Entity(:).ElectrodeID] == Analog_dataChannels(analog_channel_num));
                    Entitynull(EntityIdx2) = 1;
                end

                %             EntityIndices(analog_channel_num) = find([hFile.Entity(:).ElectrodeID] == Analog_dataChannels(analog_channel_num));
                %             fileTypeNum(analog_channel_num) = hFile.Entity(EntityIndices(analog_channel_num)).FileType;
                %             fileType{analog_channel_num} = hFile.FileInfo(1).Type;
                %             entityID(analog_channel_num) = EntityIndices(analog_channel_num);
                %             % Extract channel info
                %             [ns_RESULT, entityInfo(analog_channel_num)] = ns_GetEntityInfo(hFile, entityID(end));
            end

            %% Extract all 30k or 1k analog data channels (in basic Ripple setup, there are only 4)
            disp('Extracting Analog Data')
            if(tail_idx ~= 1)
                tmp = EMG_Data_Time_s;
                tmp2 = EMG_Data_mV;
            end
            EMG_Data_Time_s = [];
            EMG_Data_mV = [];
            EMG_Data_Time_s2 = [];
            EMG_Data_mV2 = [];

            emginfo.scale(find(Entitynull == 1)) = [hFile.Entity(EntityIndices(:)).Scale];
            fent = find(Entitynull == 1);
            for entIdx = 1:length(find(Entitynull == 1))
                emginfo.units{fent(entIdx)} = [hFile.Entity(EntityIndices(entIdx)).Units];
            end

            if(tail_idx == 1)
                nTimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices(1)).FileType).TimeStamps;
                [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices(1));
                [ns_RESULT, EMG_Data_mV_tmp] = ns_GetAnalogDataBlock(hFile, EntityIndices, 1, nTimeStamps(2,1), 0);

                EMG_Data_mV_tmp2 = zeros(size(EMG_Data_mV_tmp, 1), numAnalog_Channels);
                EMG_Data_mV_tmp2(:, find(Entitynull == 1)) = EMG_Data_mV_tmp;

                EMG_Data_Time_s(1:numAnalog_Channels,:) = repmat((0:nTimeStamps(2,1)-1) ./ analogInfo.SampleRate, [numAnalog_Channels, 1]);
                EMG_Data_mV(1:numAnalog_Channels,:) = EMG_Data_mV_tmp2';
            else
                nTimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices(1)).FileType).TimeStamps;
                [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices(1));
                [ns_RESULT, EMG_Data_mV_tmp] = ns_GetAnalogDataBlock(hFile, EntityIndices, 1, nTimeStamps(2,1), 0);

                EMG_Data_mV_tmp2 = zeros(size(EMG_Data_mV_tmp, 1), numAnalog_Channels);
                EMG_Data_mV_tmp2(:, find(Entitynull == 1)) = EMG_Data_mV_tmp;

                EMG_Data_Time_s2(1:numAnalog_Channels,:) = repmat((0:nTimeStamps(2,1)-1) ./ analogInfo.SampleRate, [numAnalog_Channels, 1]);
                EMG_Data_mV2(1:numAnalog_Channels,:) = EMG_Data_mV_tmp2';

                EMG_Data_Time_s(1:length(EntityIndices),:) = [tmp(:, :), EMG_Data_Time_s2 + time_offset];
                EMG_Data_mV(1:length(EntityIndices),:) = [tmp2(:, :), EMG_Data_mV2];
            end

            % __________________________________________________________________

%             AcuteEMG_ArrayData_Stream = {};
%             if(Good_AcuteEMG_Electrodes(1))
%                 AcuteEMG_ArrayData_Stream{1,1}(1,:) = EMG_Data_mV(1,:); % EMG 1 data in mV
%                 AcuteEMG_ArrayData_Stream{1,1}(2,:) = EMG_Data_Time_s(1,:); % timestamps (depends on 1K vs. 30K)
%                 Info2.EMGChannel{1,1} = 'EMG1';
%             end
%             if(Good_AcuteEMG_Electrodes(2))
%                 AcuteEMG_ArrayData_Stream{2,1}(1,:) = EMG_Data_mV(2,:); % EMG 1 data in mV
%                 AcuteEMG_ArrayData_Stream{2,1}(2,:) = EMG_Data_Time_s(2,:); % timestamps (depends on 1K vs. 30K)
%                 Info2.EMGChannel{2,1} = 'EMG2';
%             end
%             if(Good_AcuteEMG_Electrodes(3))
%                 AcuteEMG_ArrayData_Stream{3,1}(1,:) = EMG_Data_mV(3,:); % EMG 1 data in mV
%                 AcuteEMG_ArrayData_Stream{3,1}(2,:) = EMG_Data_Time_s(3,:); % timestamps (depends on 1K vs. 30K)
%                 Info2.EMGChannel{3,1} = 'EMG3';
%             end

        else
            AcuteEMG_ArrayData_Stream = {};
        end

        if(tail_idx == length(tail_list))
            AnalogData.data = Analog_Data_mV;
            if(Acute_EMG_DATA); AnalogData.emgdata = EMG_Data_mV; end
            AnalogData.emginfo.SampleRate = 30000; 
%             AnalogData.PD_ArrayData_Stream = PD_ArrayData_Stream;
%             AnalogData.AcuteEMG_ArrayData_Stream = AcuteEMG_ArrayData_Stream;
%             AnalogData.AcuteEMG_info = emginfo;
            analoginfo.SampleRate = Info.SampleRate;
            analoginfo.Channels = {'EyeX', 'EyeY', 'Pupil', 'PhotoDiode Aligned to First Peripheral Target Onset'};
            
            AnalogData.AnalogInfo = analoginfo;
        end
    else
        disp('Did NOT get Analog Eye Channels for Behavioural Analysis')
        disp('Did NOT get Analog Photodiode Channel(s) for Behavioural Analysis')
    end

    %% Analog Potted Link R32 Data you wish to extract
    if R32_EMG_DATA >= 1
        disp('Getting Potted Link Data for EMG Analysis')

        % General Recording Info
        numEMG_Channels = 32; % There are 32 channels (not including ground/reference) that come from the R32

        % Complete Data File Path
        completeFilePath_EMGAnalog = strcat(Folders.O_NHP_Folder,DataPath,DataFolder,'/',DataFileName,'.ns3');

        %% Open the analog file and extract some basic information
        disp('opening Chronic EMG datafile')
        [ns_status, hFile] = ns_OpenFile(completeFilePath_EMGAnalog,'single'); % opens only the .ns3 file
        if strcmp(ns_status,'ns_OK')
            disp('Chronic EMG Datafile opened properly')
        else
            disp('Did not open EMG datafile properly')
        end
        %% Get the EnitityIndicies ("EntityIDs") for the EMG Channels
        for channel_num = 1:numEMG_Channels
            datafile_channel_num = channel_num + 128; % see [ns_status, hFile] = ns_OpenFile() as to why this works
            EntityIndices(channel_num) = find([hFile.Entity(:).ElectrodeID] == datafile_channel_num);
            [ns_RESULT, entityInfo(channel_num)] = ns_GetEntityInfo(hFile, EntityIndices(channel_num));
        end

        %% Extract all EMG data channels
        disp('Extracting Chronic EMG Data')
        for channel_num_AnalogEMGData = 1:length(EntityIndices)
            [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices(channel_num_AnalogEMGData));     % analog info contains things like range and sampling rate

            TimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices(channel_num_AnalogEMGData)).FileType).TimeStamps;
            numSamples = sum(TimeStamps(:,end));
            analogInputData = zeros(1,numSamples);
            startIndex = 1;
            indexCount = TimeStamps(2,1);
            for i = 1:size(TimeStamps,2)
                [~, ~, tempData] = ns_GetAnalogData(hFile, EntityIndices(channel_num_AnalogEMGData), startIndex, indexCount);
                dataRange = TimeStamps(1,i) + (1:TimeStamps(2,i));
                analogInputData(dataRange) = tempData';
                clear tempData
                if i ~= size(TimeStamps,2)
                    startIndex = startIndex + TimeStamps(2,i);
                    indexCount = TimeStamps(2,i+1);
                end
            end

            analogInputDataTime_s = (0:numSamples-1)' ./ analogInfo.SampleRate;
            analogInputData_mV = analogInputData; clear analogInputData;
            EMG_Data_Time_s(channel_num_AnalogEMGData,:) = analogInputDataTime_s';
            EMG_Data_mV(channel_num_AnalogEMGData,:) = analogInputData_mV;
        end

        %% Convert EMG_Data_mV (millivolts) to EMG_Array_Data_mV (microvolts)
        % NOTE: multiplier converts to microVolts
        for nCh = 1:numEMG_Channels
            EMG_Array_Data_mV{nCh,1}(1,:) = (EMG_Data_mV(nCh,:))*units_multiplier;
            EMG_Array_Time{nCh,1}(2,:) = EMG_Data_Time_s(nCh,:);
        end

        %% Differentials between channels
        diff_muscle_Array_raw(1,:) = EMG_Array_Data_mV(1,:)-EMG_Array_Data_mV(2,:);
        diff_muscle_Array_raw(2,:) = EMG_Array_Data_mV(3,:)-EMG_Array_Data_mV(4,:);
        diff_muscle_Array_raw(3,:) = EMG_Array_Data_mV(5,:)-EMG_Array_Data_mV(6,:);
        diff_muscle_Array_raw(4,:) = EMG_Array_Data_mV(7,:)-EMG_Array_Data_mV(8,:);
        diff_muscle_Array_raw(5,:) = EMG_Array_Data_mV(9,:)-EMG_Array_Data_mV(10,:);
        diff_muscle_Array_raw(6,:) = EMG_Array_Data_mV(11,:)-EMG_Array_Data_mV(12,:);
        diff_muscle_Array_raw(7,:) = EMG_Array_Data_mV(13,:)-EMG_Array_Data_mV(14,:);
        diff_muscle_Array_raw(8,:) = EMG_Array_Data_mV(15,:)-EMG_Array_Data_mV(16,:);
        diff_muscle_Array_raw(9,:) = EMG_Array_Data_mV(17,:)-EMG_Array_Data_mV(18,:);
        diff_muscle_Array_raw(10,:) = EMG_Array_Data_mV(19,:)-EMG_Array_Data_mV(20,:);
        diff_muscle_Array_raw(11,:) = EMG_Array_Data_mV(21,:)-EMG_Array_Data_mV(22,:);
        diff_muscle_Array_raw(12,:) = EMG_Array_Data_mV(23,:)-EMG_Array_Data_mV(24,:);
        diff_muscle_Array_raw(13,:) = EMG_Array_Data_mV(25,:)-EMG_Array_Data_mV(26,:);
        diff_muscle_Array_raw(14,:) = EMG_Array_Data_mV(27,:)-EMG_Array_Data_mV(28,:);
        diff_muscle_Array_raw(15,:) = EMG_Array_Data_mV(29,:)-EMG_Array_Data_mV(30,:);
        diff_muscle_Array_raw(16,:) = EMG_Array_Data_mV(31,:)-EMG_Array_Data_mV(32,:);

        %% Get other R32 Potted Link Information
        disp('Extracting Other Potted Link Data')
        %% Get the EnitityIndicies ("EntityIDs" for the Other Potted Link Info
        for channel_num_other = 1:15
            datafile_channel_num_other = channel_num_other+20609; % see [ns_status, hFile] = ns_OpenFile() as to why this works
            EntityIndices_other(channel_num_other) = find([hFile.Entity(:).ElectrodeID] == datafile_channel_num_other);
            [ns_RESULT, entityInfo(channel_num_other)] = ns_GetEntityInfo(hFile, EntityIndices_other(channel_num_other));
        end

        %% Extract Other Potted Link data channels
        for channel_num_AnalogOtherData = 1:length(EntityIndices_other)
            [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices_other(channel_num_AnalogOtherData));     % analog info contains things like range and sampling rate

            TimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices_other(channel_num_AnalogOtherData)).FileType).TimeStamps;
            numSamples = sum(TimeStamps(:,end));
            analogInputData = zeros(1,numSamples);
            startIndex = 1;
            indexCount = TimeStamps(2,1);
            for i = 1:size(TimeStamps,2)
                [~, ~, tempData] = ns_GetAnalogData(hFile, EntityIndices_other(channel_num_AnalogOtherData), startIndex, indexCount);
                dataRange = TimeStamps(1,i) + (1:TimeStamps(2,i));
                analogInputData(dataRange) = tempData';
                clear tempData
                if i ~= size(TimeStamps,2)
                    startIndex = startIndex + TimeStamps(2,i);
                    indexCount = TimeStamps(2,i+1);
                end
            end

            analogInputDataTime_s = (0:numSamples-1)' ./ analogInfo.SampleRate;
            analogInputData_mV = analogInputData; clear analogInputData;
            OtherR32_Array_Data_Time_s(channel_num_AnalogOtherData,:) = analogInputDataTime_s';
            OtherR32_Array_Data(channel_num_AnalogOtherData,:) = analogInputData_mV;
        end

        for nCh_other = 1:channel_num_other
            OtherR32_Data{nCh_other,1}(1,:) = OtherR32_Array_Data(nCh_other,:);
            OtherR32_Time{nCh_other,1}(2,:) = OtherR32_Array_Data_Time_s(nCh_other,:);
        end

        %% Information Regarding Muscle Channels and Sampling Rate

        Info.Monopolar_EMG{1,1} = 'CH_1';Info.Monopolar_EMG{1,2} = 'Unknown';
        Info.Monopolar_EMG{2,1} = 'CH_2';Info.Monopolar_EMG{2,2} = 'Unknown';
        Info.Monopolar_EMG{3,1} = 'CH_3';Info.Monopolar_EMG{3,2} = 'Unknown';
        Info.Monopolar_EMG{4,1} = 'CH_4';Info.Monopolar_EMG{4,2} = 'Unknown';
        Info.Monopolar_EMG{5,1} = 'CH_5';Info.Monopolar_EMG{5,2} = 'Unknown';
        Info.Monopolar_EMG{6,1} = 'CH_6';Info.Monopolar_EMG{6,2} = 'Unknown';
        Info.Monopolar_EMG{7,1} = 'CH_7';Info.Monopolar_EMG{7,2} = 'Unknown';
        Info.Monopolar_EMG{8,1} = 'CH_8';Info.Monopolar_EMG{8,2} = 'Unknown';
        Info.Monopolar_EMG{9,1} = 'CH_9';Info.Monopolar_EMG{9,2} = 'Unknown';
        Info.Monopolar_EMG{10,1} = 'CH_10';Info.Monopolar_EMG{10,2} = 'Unknown';
        Info.Monopolar_EMG{11,1} = 'CH_11';Info.Monopolar_EMG{11,2} = 'Unknown';
        Info.Monopolar_EMG{12,1} = 'CH_12';Info.Monopolar_EMG{12,2} = 'Unknown';
        Info.Monopolar_EMG{13,1} = 'CH_13';Info.Monopolar_EMG{13,2} = 'Unknown';
        Info.Monopolar_EMG{14,1} = 'CH_14';Info.Monopolar_EMG{14,2} = 'Unknown';
        Info.Monopolar_EMG{15,1} = 'CH_15';Info.Monopolar_EMG{15,2} = 'Unknown';
        Info.Monopolar_EMG{16,1} = 'CH_16';Info.Monopolar_EMG{16,2} = 'Unknown';
        Info.Monopolar_EMG{17,1} = 'CH_17';Info.Monopolar_EMG{17,2} = 'Unknown';
        Info.Monopolar_EMG{18,1} = 'CH_18';Info.Monopolar_EMG{18,2} = 'Unknown';
        Info.Monopolar_EMG{19,1} = 'CH_19';Info.Monopolar_EMG{19,2} = 'Unknown';
        Info.Monopolar_EMG{20,1} = 'CH_20';Info.Monopolar_EMG{20,2} = 'Unknown';
        Info.Monopolar_EMG{21,1} = 'CH_21';Info.Monopolar_EMG{21,2} = 'Unknown';
        Info.Monopolar_EMG{22,1} = 'CH_22';Info.Monopolar_EMG{22,2} = 'Unknown';
        Info.Monopolar_EMG{23,1} = 'CH_23';Info.Monopolar_EMG{23,2} = 'Unknown';
        Info.Monopolar_EMG{24,1} = 'CH_24';Info.Monopolar_EMG{24,2} = 'Unknown';
        Info.Monopolar_EMG{25,1} = 'CH_25';Info.Monopolar_EMG{25,2} = 'Unknown';
        Info.Monopolar_EMG{26,1} = 'CH_26';Info.Monopolar_EMG{26,2} = 'Unknown';
        Info.Monopolar_EMG{27,1} = 'CH_27';Info.Monopolar_EMG{27,2} = 'Unknown';
        Info.Monopolar_EMG{28,1} = 'CH_28';Info.Monopolar_EMG{28,2} = 'Unknown';
        Info.Monopolar_EMG{29,1} = 'CH_29';Info.Monopolar_EMG{29,2} = 'Unknown';
        Info.Monopolar_EMG{30,1} = 'CH_30';Info.Monopolar_EMG{30,2} = 'Unknown';
        Info.Monopolar_EMG{31,1} = 'CH_31';Info.Monopolar_EMG{31,2} = 'Unknown';
        Info.Monopolar_EMG{32,1} = 'CH_32';Info.Monopolar_EMG{32,2} = 'Unknown';

        Info.Bipolar_EMG{1,1} = 'CH_1';Info.Bipolar_EMG{1,2} = 'Unknown';
        Info.Bipolar_EMG{2,1} = 'CH_2';Info.Bipolar_EMG{2,2} = 'Unknown';
        Info.Bipolar_EMG{3,1} = 'CH_3';Info.Bipolar_EMG{3,2} = 'Unknown';
        Info.Bipolar_EMG{4,1} = 'CH_4';Info.Bipolar_EMG{4,2} = 'Unknown';
        Info.Bipolar_EMG{5,1} = 'CH_5';Info.Bipolar_EMG{5,2} = 'Unknown';
        Info.Bipolar_EMG{6,1} = 'CH_6';Info.Bipolar_EMG{6,2} = 'Unknown';
        Info.Bipolar_EMG{7,1} = 'CH_7';Info.Bipolar_EMG{7,2} = 'Unknown';
        Info.Bipolar_EMG{8,1} = 'CH_8';Info.Bipolar_EMG{8,2} = 'Unknown';
        Info.Bipolar_EMG{9,1} = 'CH_9';Info.Bipolar_EMG{9,2} = 'Unknown';
        Info.Bipolar_EMG{10,1} = 'CH_10';Info.Bipolar_EMG{10,2} = 'Unknown';
        Info.Bipolar_EMG{11,1} = 'CH_11';Info.Bipolar_EMG{11,2} = 'Unknown';
        Info.Bipolar_EMG{12,1} = 'CH_12';Info.Bipolar_EMG{12,2} = 'Unknown';
        Info.Bipolar_EMG{13,1} = 'CH_13';Info.Bipolar_EMG{13,2} = 'Unknown';
        Info.Bipolar_EMG{14,1} = 'CH_14';Info.Bipolar_EMG{14,2} = 'Unknown';
        Info.Bipolar_EMG{15,1} = 'CH_15';Info.Bipolar_EMG{15,2} = 'Unknown';
        Info.Bipolar_EMG{16,1} = 'CH_16';Info.Bipolar_EMG{16,2} = 'Unknown';

        Info.SamplingRate = 2000;

        ChronicEMGData.Monopolar_Raw = EMG_Array_Data_mV;
        ChronicEMGData.EMG_Time = EMG_Array_Time;
        ChronicEMGData.Bipolar_Raw = diff_muscle_Array_raw;
        ChronicEMGData.OtherR32_Time = OtherR32_Time;
        ChronicEMGData.OtherR32_Data = OtherR32_Data;
        ChronicEMGData.Info = Info;
    else
        disp('Did NOT get Potted Link Data for EMG Analysis')
    end


    %% Extract Neural Data
    if Neural_DATA == 1
        if FHC_DATA == 1
            disp('You have NOT written the extraction code for FHC extraction yet!')
        else
            disp('Did not extract FHC data')
        end
        if Neural_Array1_Data == 1
            disp('You have NOT written the extraction code for Neural_Array1_Data yet!')
        else
            disp('Did not extract Neural_Array1_Data data')
        end
        if Neural_Array2_Data == 1
            disp('You have NOT written the extraction code for Neural_Array2_Data yet!')
        else
            disp('Did not extract Neural_Array2_Data data')
        end

        if Plexon_Data == 1
            disp('Extracting 16 (or 32) Channels of Raw Data from Plexon Electrode')
            % General Recording Info
            numNeural_Channels = numberOfChan; % There are 32 possible channels that come from the Plexon S-Probe

            % Complete Data File Path
            completeFilePath_Plexon = strcat(Folders.O_NHP_Folder,DataPath,DataFolder,'/',DataFileName,'.ns5');

            %% Open the analog file and extract some basic information
            disp('opening Neural datafile')
            [ns_status, hFile] = ns_OpenFile(completeFilePath_Plexon,'single'); % opens only the .ns5 file
            if strcmp(ns_status,'ns_OK')
                disp('Raw Neural Datafile opened properly')
            else
                disp('Did not open Raw Neural datafile properly')
            end

            %% Get the EnitityIndicies ("EntityIDs") for the Neural Channels

            for channel_num = (ch_offset+1):(ch_offset+numNeural_Channels)
                datafile_channel_num = channel_num; % see [ns_status, hFile] = ns_OpenFile() as to why this works.
                
                EntityIndices(channel_num - ch_offset) = find(([hFile.Entity(:).ElectrodeID]) == datafile_channel_num);
                [ns_RESULT, entityInfo(channel_num)] = ns_GetEntityInfo(hFile, EntityIndices(channel_num - ch_offset));
            end


            %             for channel_num = 1:numNeural_Channels
            %                 %         channel_num
            %                 %         hFile.Entity(:).ElectrodeID
            %
            %                 datafile_channel_num = channel_num; % see [ns_status, hFile] = ns_OpenFile() as to why this works
            %                 EntityIndices(channel_num) = find(([hFile.Entity(:).ElectrodeID] - ch_offset) == datafile_channel_num);
            %                 [ns_RESULT, entityInfo(channel_num)] = ns_GetEntityInfo(hFile, EntityIndices(channel_num));
            %             end

            %% Extract all raw Neural data channels
            disp('Extracting Raw Neural Data')
            if(tail_idx ~= 1)
                tmp = RawNeural_Data_Time_s;
                tmp2 = RawNeural_Data_mV;
            end
            RawNeural_Data_Time_s = [];
            RawNeural_Data_mV = [];
            RawNeural_Data_Time_s2 = [];
            RawNeural_Data_mV2 = [];

            rawinfo.scale = [hFile.Entity(EntityIndices(1:32)).Scale];
            for entIdx = 1:32
                rawinfo.units{entIdx} = [hFile.Entity(EntityIndices(entIdx)).Units];
            end

            if(tail_idx == 1)
                nTimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices(1)).FileType).TimeStamps;
                [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices(1));
                [ns_RESULT, RawNeural_Data_mV_tmp] = ns_GetAnalogDataBlock(hFile, EntityIndices, 1, nTimeStamps(2,1), 0);
                RawNeural_Data_Time_s(1:length(EntityIndices),:) = repmat((0:nTimeStamps(2,1)-1) ./ analogInfo.SampleRate, [length(EntityIndices), 1]);
                RawNeural_Data_mV(1:length(EntityIndices),:) = RawNeural_Data_mV_tmp';
            else
                nTimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices(1)).FileType).TimeStamps;
                [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices(1));
                [ns_RESULT, RawNeural_Data_mV_tmp] = ns_GetAnalogDataBlock(hFile, EntityIndices, 1, nTimeStamps(2,1), 0);
                RawNeural_Data_Time_s2(1:length(EntityIndices),:) = repmat((0:nTimeStamps(2,1)-1) ./ analogInfo.SampleRate, [length(EntityIndices), 1]);
                RawNeural_Data_mV2(1:length(EntityIndices),:) = RawNeural_Data_mV_tmp';
                RawNeural_Data_Time_s(1:length(EntityIndices),:) = [tmp(:, :), RawNeural_Data_Time_s2 + time_offset];
                RawNeural_Data_mV(1:length(EntityIndices),:) = [tmp2(:, :), RawNeural_Data_mV2];
            end

            rawinfo.SampleRate = analogInfo.SampleRate;
            %             for channel_num_RawNeuralData = 1:length(EntityIndices)
            %
            %                 electrode_number = hFile.Entity(EntityIndices(channel_num_RawNeuralData)).ElectrodeID - ch_offset;
            %                 units = hFile.Entity(EntityIndices(channel_num_RawNeuralData)).Units;
            %                 scale = hFile.Entity(EntityIndices(channel_num_RawNeuralData)).Scale;
            %                 label = hFile.Entity(EntityIndices(channel_num_RawNeuralData)).Label;
            %                 count = hFile.Entity(EntityIndices(channel_num_RawNeuralData)).Count;
            %
            %                 [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices(channel_num_RawNeuralData));     % analog info contains things like range and sampling rate
            %
            %                 TimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices(channel_num_RawNeuralData)).FileType).TimeStamps;
            %                 numSamples = sum(TimeStamps(:,end));
            %                 analogInputData = zeros(1,numSamples);
            %                 startIndex = 1;
            %                 indexCount = TimeStamps(2,1);
            %                 for i = 1:size(TimeStamps,2)
            %                     [~, ~, tempData] = ns_GetAnalogData(hFile, EntityIndices(channel_num_RawNeuralData), startIndex, indexCount);
            %                     dataRange = TimeStamps(1,i) + (1:TimeStamps(2,i));
            %                     analogInputData(dataRange) = tempData';
            %                     if i ~= size(TimeStamps,2)
            %                         startIndex = startIndex + TimeStamps(2,i);
            %                         indexCount = TimeStamps(2,i+1);
            %                     end
            %                 end
            %
            %                 analogInputDataTime_s = (0:numSamples-1)' ./ analogInfo.SampleRate;
            %                 analogInputData_mV = analogInputData;
            %                 RawNeural_Data_Time_s(channel_num_RawNeuralData,:) = analogInputDataTime_s';
            %                 RawNeural_Data_mV(channel_num_RawNeuralData,:) = analogInputData_mV;
            %
            %                 RAWNEURALData.Info(channel_num_RawNeuralData).ElectrodeNumber = electrode_number;
            %                 RAWNEURALData.Info(channel_num_RawNeuralData).Units = units;
            %                 RAWNEURALData.Info(channel_num_RawNeuralData).Scale = scale;
            %                 RAWNEURALData.Info(channel_num_RawNeuralData).Label = label;
            %                 RAWNEURALData.Info(channel_num_RawNeuralData).Count = count;
            %                 RAWNEURALData.Info(channel_num_RawNeuralData).SampleRate = 30000;
            %             end

            % Whole Data Stream Array for Saving
%             RAWNEURAL_ArrayData_Stream = {};
% 
%             for nChRawNeural = 1:numNeural_Channels
%                 RAWNEURAL_ArrayData_Stream{nChRawNeural,1}(1,:) = RawNeural_Data_mV(nChRawNeural,:);
%                 RAWNEURAL_ArrayData_Stream{nChRawNeural,1}(2,:) = RawNeural_Data_Time_s(nChRawNeural,:);
%             end
% 
%             RAWNEURALData.RAW_Data_Stream = RAWNEURAL_ArrayData_Stream;
%             RAWNEURALData.info = rawinfo;

            %% Getting Digital Timestamps for Neural Data

            % Complete Data File Path: Digital Data
            completeFilePath_Digital = strcat(Folders.O_NHP_Folder,DataPath,DataFolder,'/',[DataFileName],'.nev');%,[DataFileName, '_Allsorted'],'.nev');

            % Open the file and extract some basic information
            %     disp('opening Digital datafile to get neural time stamps')
            %     [ns_status, hFile] = ns_OpenFile(completeFilePath_Digital,'single');
            %         if strcmp(ns_status,'ns_OK')
            %             disp('Digital datafile for neural timestamps opened properly')
            %         else
            %             disp('Did not open Digital datafile for neural time stamps properly')
            %         end
            %
            %     for nNeuralEvents = 1:length(hFile.Entity)
            %         entityID = nNeuralEvents;
            %         electrode_number = hFile.Entity(nNeuralEvents).ElectrodeID;
            %         units = hFile.Entity(nNeuralEvents).Units;
            %         scale = hFile.Entity(nNeuralEvents).Scale;
            %         nUnits = hFile.Entity(nNeuralEvents).nUnits;
            %         label = hFile.Entity(nNeuralEvents).Label;
            %
            %         if electrode_number ~= 0 % this is the parallel port inputs (i.e., strobecodes)
            %         [ns_RESULT, nsSegmentInfo] = ns_GetSegmentInfo(hFile, entityID);
            %         numCount = hFile.Entity(nNeuralEvents).Count;
            %         samplerate = nsSegmentInfo.SampleRate;
            %
            %             neural_timeStamps = NaN(1,numCount);
            % %             for a = 1:numCount
            % %                 voltage_data{a} = NaN;
            % %             end
            %             SampleCount = NaN(1,numCount);
            %             UnitID = NaN(1,numCount);
            % %             numCount
            %             for i = 1:numCount
            %                 %[~, neural_timeStamps, neural_data, SampleCount, UnitID] = ns_GetSegmentData(hFile, entityID(end), i);
            %                 [~, NTS, nd, SC, UID] = ns_GetSegmentData(hFile, entityID, i);
            % %                 NTS
            %                 neural_timeStamps(1,i) = NTS;
            % %                 voltage_data{i} = nd;
            %                 SampleCount(1,i) = SC;
            %                 UnitID(1,i) = UID;
            %             end
            % %             clear ChannelID;
            %             ChannelID = electrode_number;
            %             spkData(ChannelID).neural_timeStamps = neural_timeStamps;
            % %             spkData(ChannelID).voltage_data = voltage_data;
            %             spkData(ChannelID).SampleCount = SampleCount;
            %             spkData(ChannelID).UnitID = UnitID;
            %             spkData(ChannelID).ElectrodeNumber = electrode_number;
            %             spkData(ChannelID).Units = units;
            %             spkData(ChannelID).Scale = scale;
            %             spkData(ChannelID).nUnits = nUnits;
            %             spkData(ChannelID).Label = label;
            %             spkData(ChannelID).SampleRate = samplerate;
            % %                 spkData(nNeuralEvents).neural_timeStamps = neural_timeStamps;
            % %                 spkData(nNeuralEvents).voltage_data = voltage_data;
            % %                 spkData(nNeuralEvents).SampleCount = SampleCount;
            % %                 spkData(nNeuralEvents).UnitID = UnitID;
            % %                 spkData(nNeuralEvents).ElectrodeNumber = electrode_number;
            % %                 spkData(nNeuralEvents).Units = units;
            % %                 spkData(nNeuralEvents).Scale = scale;
            % %                 spkData(nNeuralEvents).nUnits = nUnits;
            % %                 spkData(nNeuralEvents).Label = label;
            % %                 spkData(nNeuralEvents).SampleRate = samplerate;
            %         end
            %     end

            %% Get LFP DATA
            % Complete Data File Path: LFP Data
            completeFilePath_Digital = strcat(Folders.O_NHP_Folder,DataPath,DataFolder,'/',DataFileName,'.ns2');

            % Open the file and extract some basic information
            disp('opening analog datafile to get LFPs')
            [ns_status, hFile] = ns_OpenFile(completeFilePath_Digital,'single');
            if strcmp(ns_status,'ns_OK')
                disp('Analog datafile opened properly to get LFPs')
            else
                disp('Did not open analog datafile properly to get LFPs')
            end

            numLFP_Channels = numberOfChan; % possible number of pins in the plexon omnetics connector for the 16ch v-probe
            % NOTE: Channels that have viable data will match up with RIPPLE channel
            % mapping file that was used during data collection!!!
            %% Get the EnitityIndicies ("EntityIDs") for the EMG Channels
            %             for channel_num = 1:numLFP_Channels
            %                 datafile_channel_num = channel_num; % see [ns_status, hFile] = ns_OpenFile() as to why this works
            % %                 [hFile.Entity(:).ElectrodeID]
            %                 EntityIndices(channel_num) = find(([hFile.Entity(:).ElectrodeID] - ch_offset) == datafile_channel_num);
            %                 [ns_RESULT, entityInfo(channel_num)] = ns_GetEntityInfo(hFile, EntityIndices(channel_num));
            %             end

            for channel_num = (ch_offset+1):(ch_offset+numNeural_Channels)
                datafile_channel_num = channel_num; % see [ns_status, hFile] = ns_OpenFile() as to why this works.
                EntityIndices(channel_num - ch_offset) = find(([hFile.Entity(:).ElectrodeID]) == datafile_channel_num);
                [ns_RESULT, entityInfo(channel_num)] = ns_GetEntityInfo(hFile, EntityIndices(channel_num - ch_offset));
            end

            %% Extract all LFP data channels
            disp('Extracting LFP Data')
            if(tail_idx ~= 1)
                tmp = LFP_Data_Time_s;
                tmp2 = LFP_Data_mV;
            end

            LFP_Data_Time_s = [];
            LFP_Data_mV = [];
            LFP_Data_Time_s2 = [];
            LFP_Data_mV2 = [];

            lfpinfo.scale = [hFile.Entity(EntityIndices(1:32)).Scale];
            for entIdx = 1:32
                lfpinfo.units{entIdx} = [hFile.Entity(EntityIndices(entIdx)).Units];
            end

            if(tail_idx == 1)
                nTimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices(1)).FileType).TimeStamps;
                [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices(1));
                [ns_RESULT, LFP_Data_mV_tmp] = ns_GetAnalogDataBlock(hFile, EntityIndices, 1, nTimeStamps(2,1), 0);
                LFP_Data_Time_s(1:length(EntityIndices),:) = repmat((0:nTimeStamps(2,1)-1) ./ analogInfo.SampleRate, [length(EntityIndices), 1]);
                LFP_Data_mV(1:length(EntityIndices),:) = LFP_Data_mV_tmp';
            else
                nTimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices(1)).FileType).TimeStamps;
                [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices(1));
                [ns_RESULT, LFP_Data_mV_tmp] = ns_GetAnalogDataBlock(hFile, EntityIndices, 1, nTimeStamps(2,1), 0);
                LFP_Data_Time_s2(1:length(EntityIndices),:) = repmat((0:nTimeStamps(2,1)-1) ./ analogInfo.SampleRate, [length(EntityIndices), 1]);
                LFP_Data_mV2(1:length(EntityIndices),:) = LFP_Data_mV_tmp';
                LFP_Data_Time_s(1:length(EntityIndices),:) = [tmp(:, :), LFP_Data_Time_s2 + time_offset];
                LFP_Data_mV(1:length(EntityIndices),:) = [tmp2(:, :), LFP_Data_mV2];
            end
            
            lfpinfo.SampleRate = analogInfo.SampleRate;

            %             for channel_num_LFPData = 1:length(EntityIndices)
            %
            %                 %             electrode_number = hFile.Entity(EntityIndices(channel_num_LFPData)).ElectrodeID;
            %                 %             units = hFile.Entity(EntityIndices(channel_num_LFPData)).Units;
            %                 %             scale = hFile.Entity(EntityIndices(channel_num_LFPData)).Scale;
            %                 %             label = hFile.Entity(EntityIndices(channel_num_LFPData)).Label;
            %                 %             count = hFile.Entity(EntityIndices(channel_num_LFPData)).Count;
            %
            %                 [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, EntityIndices(channel_num_LFPData));     % analog info contains things like range and sampling rate
            %                 TimeStamps = hFile.FileInfo(hFile.Entity(EntityIndices(channel_num_LFPData)).FileType).TimeStamps;
            %                 numSamples = sum(TimeStamps(:,end));
            %                 analogInputData = zeros(1,numSamples);
            %                 startIndex = 1;
            %                 indexCount = TimeStamps(2,1);
            %                 for i = 1:size(TimeStamps,2)
            %                     [~, ~, tempData] = ns_GetAnalogData(hFile, EntityIndices(channel_num_LFPData), startIndex, indexCount);
            %                     dataRange = TimeStamps(1,i) + (1:TimeStamps(2,i));
            %                     analogInputData(dataRange) = tempData';
            %                     if i ~= size(TimeStamps,2)
            %                         startIndex = startIndex + TimeStamps(2,i);
            %                         indexCount = TimeStamps(2,i+1);
            %                     end
            %                 end
            %                 analogInputDataTime_s = (0:numSamples-1)' ./ analogInfo.SampleRate;
            %                 analogInputData_mV = analogInputData;
            %                 LFP_Data_Time_s(channel_num_LFPData,:) = analogInputDataTime_s';
            %                 LFP_Data_mV(channel_num_LFPData,:) = analogInputData_mV;
            %
            %                 %             LFPData.Info(channel_num_LFPData).ElectrodeNumber = electrode_number;
            %                 %             LFPData.Info(channel_num_LFPData).Units = units;
            %                 %             LFPData.Info(channel_num_LFPData).Scale = scale;
            %                 %             LFPData.Info(channel_num_LFPData).Label = label;
            %                 %             LFPData.Info(channel_num_LFPData).Count = count;
            %                 LFPData.Info(channel_num_LFPData).SampleRate = 1000;
            %             end

            %% Whole LFP Data Stream Array

%             LFP_ArrayData_Stream = {};
% 
%             for nChLFP = 1: numLFP_Channels
%                 LFP_ArrayData_Stream{nChLFP,1}(1,:) = LFP_Data_mV(nChLFP,:);
%                 LFP_ArrayData_Stream{nChLFP,1}(2,:) = LFP_Data_Time_s(nChLFP,:);
%             end
% 
%             LFPData.LFP_Data_Stream = LFP_ArrayData_Stream;
%             LFPData.info = lfpinfo;

            disp('LFP Data Extracted Properly')

            %% Create Plexon S-Probe Structure

            if(tail_idx == length(tail_list))
                Plexon_SProbe.RAW_Data.data = RawNeural_Data_mV([1:2:32, 2:2:32], :);
                Plexon_SProbe.RAW_Data.info = rawinfo;
                %     Plexon_SProbe.Spk_Data = spkData;
                Plexon_SProbe.LFP_Data.data = LFP_Data_mV([1:2:32, 2:2:32], :);

                if(Physio_Info.numCh == 16)
                    Plexon_SProbe.LFP_Data.data = Plexon_SProbe.LFP_Data.data([1:16], :);
                    Plexon_SProbe.RAW_Data.data = Plexon_SProbe.RAW_Data.data([1:16], :);
                end

                Plexon_SProbe.LFP_Data.info = lfpinfo;
            end

        else
            disp('Did NOT get Neural Data for 16 Channel Plexon Electrode')
        end
    end

    time_offset = RawNeural_Data_Time_s(1, end);
end

%% OUTPUT ALL RAW PHYSIOLOGY DATA AS A STRUCTURE
PHYSIOLOGY_RAW =[];
% Strobe Codes and associated timestamps
if exist('DigitalInfo') == 1
    PHYSIOLOGY_RAW.DigitalInfo = DigitalInfo;
else
    PHYSIOLOGY_RAW.DigitalInfo = [];
end
% EyeX, EyeY, Pupil, and Acute EMG data
if exist('AnalogData') == 1
    PHYSIOLOGY_RAW.AnalogData = AnalogData;
else
    PHYSIOLOGY_RAW.AnalogData = [];
end
% Potted Link (Chronic) EMG
if exist('ChronicEMGData') == 1
    PHYSIOLOGY_RAW.ChronicEMGData = ChronicEMGData;
else
    PHYSIOLOGY_RAW.ChronicEMGData = [];
end
% Plexon_SProbe
if exist('Plexon_SProbe')
    % Change channels order to the correct order
    PHYSIOLOGY_RAW.Plexon_SProbe = Plexon_SProbe;
else
    PHYSIOLOGY_RAW.Plexon_SProbe = [];
end
% FHC Single Contact
if exist('FHC')
    PHYSIOLOGY_RAW.FHC = FHC;
else
    PHYSIOLOGY_RAW.FHC = [];
end
% Neural_Array1
if exist('Neural_Array1')
    PHYSIOLOGY_RAW.Neural_Array1 = Neural_Array1;
else
    PHYSIOLOGY_RAW.Neural_Array1 = [];
end
% Neural_Array2
if exist('Neural_Array2')
    PHYSIOLOGY_RAW.Neural_Array2 = Neural_Array2;
else
    PHYSIOLOGY_RAW.Neural_Array2 = [];
end

Z.PHYSIOLOGY_RAW = PHYSIOLOGY_RAW;

Z_Raw.data = Z.PHYSIOLOGY_RAW.Plexon_SProbe.RAW_Data.data;
Z_Raw.info = Z.PHYSIOLOGY_RAW.Plexon_SProbe.RAW_Data.info;

Z_LFP.data = Z.PHYSIOLOGY_RAW.Plexon_SProbe.LFP_Data.data;
Z_LFP.info = Z.PHYSIOLOGY_RAW.Plexon_SProbe.LFP_Data.info;

% What are the relevant event codes for splicing data in "Ripple_Splice_Trials"?
ST.statecode_TrialStart = 9; % Start of a trial (MonkeyLogic default)
ST.statecode_TrialEnd = 18; % End of a trial (MonkeyLogic default)
ST.statecode_Reward = 90; % Rewarded Trials (End User defined in in MonkeyLogic)
ST.statecode_CueToGo_ET_ST = 60; % Cue to Go for Emerging Target,Single Target Task (121)
ST.statecode_CueToGo_ET_DT = 200; % Cue to Go for Emerging Target, Double Target Task (2XX)
ST.statecode_CueToGo_EyeHand = 30;% Cue to Go for eye only and eye-hand  (504,506)

Z_Digital = [];
if(Physio_Info.Statecodes)
    clear Z_Digital
    Z_Digital.data = [Z.PHYSIOLOGY_RAW.DigitalInfo.Statecodes; Z.PHYSIOLOGY_RAW.DigitalInfo.Statecode_timeStamps];
    Z_Digital.info = ST;
end

Z_Analog.data = Z.PHYSIOLOGY_RAW.AnalogData.data;
Z_Analog.info = Z.PHYSIOLOGY_RAW.AnalogData.AnalogInfo;

Z_EMG = [];
if(Physio_Info.Acute_EMG_DATA)
    clear Z_EMG
    Z_EMG.data = Z.PHYSIOLOGY_RAW.AnalogData.emgdata;
    Z_EMG.info = Z.PHYSIOLOGY_RAW.AnalogData.emginfo;
end

cd(Folders.save_dir)
disp('Saving Z files ... it usually takes a while (3-7 mins) if you want to save the raw file, otherwise 20 secs! ...')

if(saveRaw); save(['Z_', ses_name_main, '_', region, '_Raw'], 'Z_Raw','-v7.3'); end
save(['Z_', ses_name_main, '_', region, '_LFP'], 'Z_LFP','-v7.3')
if(Physio_Info.Statecodes)
    save(['Z_', ses_name_main, '_Digital'], 'Z_Digital','-v7.3')
    load(['Z_', ses_name_main, '_ML'])
    PD_onset_tr = get_PD_trial(Z_ML, Z_Digital, Z_Analog);
    Z_Analog.PD_sec = PD_onset_tr;
end
save(['Z_', ses_name_main, '_Analog'], 'Z_Analog','-v7.3')
if(Physio_Info.Acute_EMG_DATA); save(['Z_', ses_name_main, '_EMG'], 'Z_EMG','-v7.3'); end





% Spikes from your sorter
% if(1) % strcmp(Physio_Info.Sorter, 'Kilosort')
%     disp('Getting Sorter data ...')
%     cd(Folders.datFolder)
%     load([ses_name_main, '_Phy_', region])
%     Z.NumberofChannels = Physio_Info.numCh;
% 
%     for i = 1:Z.NumberofChannels
%         spk_arr = [];
%         ID_arr = [];
%         ch_cluster = find(cluster_channel == i);
%         for k = 1:length(ch_cluster)
%             spk_arr = [spk_arr, spikes_phy(2, find(spikes_phy(1, :) == (clusterId(1, ch_cluster(k)))))/30000];
%             ID_arr = [ID_arr, spikes_phy(1, find(spikes_phy(1, :) == (clusterId(1, ch_cluster(k)))))];
%         end
% 
%         [spk_arr, I] = sort(spk_arr);
%         Z.PHYSIOLOGY_RAW.Plexon_SProbe.Spk_Data(i).SampleRate = 30000;
%         Z.PHYSIOLOGY_RAW.Plexon_SProbe.Spk_Data(i).neural_timeStamps = spk_arr;
%         Z.PHYSIOLOGY_RAW.Plexon_SProbe.Spk_Data(i).UnitID = ID_arr(I);
%         Z.PHYSIOLOGY_RAW.Plexon_SProbe.Spk_Data(i).nUnits = length(ch_cluster);
%         Z.num_unit_offline(i) = length(ch_cluster);
%         Z.units{i} = clusterId(1, ch_cluster);
%     end
%     Z.goodUnits = clusterId(find(clusterScore == 2));
%     Z.unitScores = clusterScore;
% 
% elseif(strcmp(Physio_Info.Sorter, 'Plexon'))
% 
%     load(['Gr', Physio_Info.DataFileName(2:end), '_Plexon_Sorted_plx'])
% 
%     for i = 1:32
% 
%         Z.PHYSIOLOGY_RAW.Plexon_SProbe.Spk_Data(i).neural_timeStamps = double(Spk_Data.neural_timeStamps{i,1}');
%         Z.PHYSIOLOGY_RAW.Plexon_SProbe.Spk_Data(i).UnitID = double(Spk_Data.UnitID{i,1})';
%         Z.PHYSIOLOGY_RAW.Plexon_SProbe.Spk_Data(i).nUnits = Spk_Data.nUnits(i);
% 
%     end
% 
%     Z.num_unit_offline = Spk_Data.nUnits;
% end


%Ripple_Sanity_Plots(Physio_Info, PHYSIOLOGY_RAW)


