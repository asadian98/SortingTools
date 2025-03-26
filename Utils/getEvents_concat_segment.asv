function getEvents_concat_segment(ses_name_main, Z, animal_name, tail, Physio_Info, Folders)
% clc; clear;

% ses_name_main = 'Be240405';
% animal_name = 'Belle';

tail_ML = tail;

event_times_l = [];
event_times_r = [];PD_Trial = [];
statecodes_all = [];
statecode_timeStamps_all = [];
PD_ArrayData_Stream_all1 = [];
PD_ArrayData_Stream_all2 = [];
PD_ArrayData_Stream_all = [];
time_offset = 0;

for tail_idx = 1:length(tail)
    ses_name = [ses_name_main, tail{tail_idx}];
    completeFilePath_Analog = [Folders.Ripple_Raw_Data_Folder, animal_name, '\', ses_name, '.ns2'];
    completeFilePath_Analog
    [ns_status, hFile] = ns_OpenFile(completeFilePath_Analog,'single');
    if strcmp(ns_status,'ns_OK')
        disp('Analog datafile opened properly')
    else
        disp('Did not open Analog datafile properly')
    end

    AnalogSampleRate = 1000;
    Analog_dataChannels = 10244; % Photodiode
    analog_channel_num = length(Analog_dataChannels);
    [hFile.Entity(:).ElectrodeID]
    EntityIndices(analog_channel_num) = find([hFile.Entity(:).ElectrodeID] == Analog_dataChannels(analog_channel_num));
    fileTypeNum(analog_channel_num) = hFile.Entity(EntityIndices(analog_channel_num)).FileType;
    fileType{analog_channel_num} = hFile.FileInfo(1).Type;
    entityID(analog_channel_num) = EntityIndices(analog_channel_num);
    % Extract channel info
    [ns_RESULT, entityInfo(analog_channel_num)] = ns_GetEntityInfo(hFile, entityID(end));

    Analog_Data_Time_s = [];
    Analog_Data_mV = [];
    for channel_num_AnalogData = 1:length(EntityIndices)
        [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, entityID(channel_num_AnalogData));     % analog info contains things like range and sampling rate

        TimeStamps = hFile.FileInfo(hFile.Entity(entityID(channel_num_AnalogData)).FileType).TimeStamps;
        numSamples = sum(TimeStamps(:,end));
        analogInputData = zeros(1,numSamples);
        startIndex = 1;
        indexCount = TimeStamps(2,1);
        for i = 1:size(TimeStamps,2)
            [~, ~, tempData] = ns_GetAnalogData(hFile, entityID(channel_num_AnalogData), startIndex, indexCount);
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
        Analog_Data_Time_s(channel_num_AnalogData,:) = analogInputDataTime_s';
        Analog_Data_mV(channel_num_AnalogData,:) = analogInputData_mV;
    end

    PD_ArrayData_Stream = [];

    PD_ArrayData_Stream(1,:) = Analog_Data_mV(1,:); % Photodiode data in mV
    PD_ArrayData_Stream(2,:) = Analog_Data_Time_s(1,:); % timestamps (depends on 1K vs. 30K)

    PD_ArrayData_Stream_all1 = [PD_ArrayData_Stream_all1, PD_ArrayData_Stream(1,:)];
    PD_ArrayData_Stream_all2 = [PD_ArrayData_Stream_all2, PD_ArrayData_Stream(2,:) + time_offset];

    PD_ArrayData_Stream_all = [PD_ArrayData_Stream_all1; PD_ArrayData_Stream_all2];

    if 1
        disp('Getting Digital Statecodes for Behavioural Analysis')
        % Complete Data File Path: Digital Data
        completeFilePath_Digital = [Folders.Ripple_Raw_Data_Folder, animal_name, '\', ses_name, '.nev'];
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

        DigitalInfo.Statecodes = statecodes;
        DigitalInfo.Statecode_timeStamps = statecode_timeStamps;
        DigitalInfo.DataSize = dataSize;

    else
        disp('Did NOT get Digital Statecodes for Behavioural Analysis')
    end
% 
%     Strobe_temp = [];
%     statecodes = [];
%     statecode_timeStamps = [];
%     Strobe_temp(:,1) = DigitalInfo.Statecodes; % codes sent from behavioural acquisition program (MonkeyLogic in ALC data)
%     Strobe_temp(:,2) = DigitalInfo.Statecode_timeStamps; % time is in seconds
%     if Strobe_temp(1,1) == Strobe_temp(2,1) % sample rate has produced two strobe codes for each code sent
%         statecodes(:,1) = Strobe_temp(1:2:end,1);
%         statecode_timeStamps(:,1) = Strobe_temp(1:2:end,2);
%         disp('Picked first of two strobe codes sent')
%         if length(statecodes) ~= length(statecode_timeStamps)
%             disp('ERROR(?):length of statecodes and statecode_timeStamps are not equal')
%         end
%     else
%         statecodes(:,1) = Strobe_temp(:,1);
%         statecode_timeStamps(:,1) = Strobe_temp(:,2);
%         disp('Picked only strobe code sent')
%         if length(statecodes) ~= length(statecode_timeStamps)
%             disp('ERROR(?):length of statecodes and statecode_timeStamps are not equal')
%         end
%     end

    % You should comment the next lines if you splitted the ripple file
    num9 = find(statecodes == 9);
    num18 = find(statecodes == 18);

    if(num9(1) ~= 1)
        statecodes(1:num9(1)-1) = [];
        statecode_timeStamps(1:num9(1)-1) = [];
    end

    if(num18(end) ~= length(statecodes))
        statecodes((num18(end) + 1) : end) = [];
        statecode_timeStamps((num18(end) + 1) : end) = [];
    end

    % For files that have less trials in ML file
    st_idx = find(statecodes == 9);
    st_idx2 = find(statecodes == 18);
    if(length(st_idx) > Z.tr_sep(tail_idx))
        statecodes(st_idx2(Z.tr_sep(tail_idx))+1:end) = [];
    end
    assert(length(find(statecodes == 9)) == Z.tr_sep(tail_idx), 'Error: check number of trials, somethings wrong');

    statecodes_all = [statecodes_all, statecodes];
    statecode_timeStamps_all = [statecode_timeStamps_all, statecode_timeStamps + time_offset];

    time_offset = PD_ArrayData_Stream_all2(end);

end


% statecode_timeStamps
% PD_ArrayData_Stream

ST.statecode_TrialStart = 9; % Start of a trial (MonkeyLogic default)
ST.statecode_TrialEnd = 18; % End of a trial (MonkeyLogic default)
statecode_TrialStart = ST.statecode_TrialStart; % Start of a trial (MonkeyLogic default)
statecode_TrialEnd = ST.statecode_TrialEnd ; % End of a trial (MonkeyLogic default)
statecode_indicies_statecode_TrialStart = find(statecodes_all == statecode_TrialStart);

statecode_times_TrialStart=[];
for t_statecode_TrialStart = 1:length(statecode_indicies_statecode_TrialStart)
    statecode_times_TrialStart(t_statecode_TrialStart) = statecode_timeStamps_all(statecode_indicies_statecode_TrialStart(t_statecode_TrialStart));
end
statecode_indicies_statecode_TrialEnd = find(statecodes_all == statecode_TrialEnd);
statecode_times_TrialEnd =[];
for t_statecode_TrialEnd = 1:length(statecode_indicies_statecode_TrialEnd)
    statecode_times_TrialEnd(t_statecode_TrialEnd) = statecode_timeStamps_all(statecode_indicies_statecode_TrialEnd(t_statecode_TrialEnd));
end
offset = length(PD_Trial);
for nTrials = 1:min(length(statecode_times_TrialStart), length(statecode_times_TrialEnd))
    PD_indicies = find(PD_ArrayData_Stream_all(2,:)>statecode_times_TrialStart(nTrials)...
        & PD_ArrayData_Stream_all(2,:)<statecode_times_TrialEnd(nTrials));

    %     PHYSIOLOGY_SPLICED.Trial(nTrials).Photodiode{1,1} = 'PhotoDiode_Aligned to Target Onset';
    %     PD_Trial(nTrials) = PD_ArrayData_Stream(1,PD_indicies);

    PD_temp_data = PD_ArrayData_Stream_all(1, PD_indicies);
    %     PD_Trial{nTrials} = statecode_times_TrialStart(nTrials) + find(PD_temp_data >= 2000,1);
    if isempty(PD_ArrayData_Stream_all(2, PD_indicies(find(PD_temp_data >= 2000,1))))
        PD_Trial(nTrials+offset) = NaN;
    else
        PD_Trial(nTrials+offset) = PD_ArrayData_Stream_all(2, PD_indicies(find(PD_temp_data >= 2000,1)));
    end
end

predur_min = 200;
predur_max = 500;
postdur = 500;

% Z = Concat_Z_Files(ses_name_main, Physio_Info, tail_ML, tail, animal_name, Folders);

Z.EyeRT = Z.EyeRT(:, 1);
Z.EyeRT(find(Z.EyeRT <= 0 | Z.EyeRT >= 250)) = 150;
Z.ArmRT(find(Z.ArmRT <= 0)) = 250;

lowRT_cutoff = 70;
highRT_cutoff = 300;


assert(length(PD_Trial) == size(Z.condition, 2), 'Error: check number of trials, somethings wrong');

conds = unique(Z.condition);
event_times_l = {};
event_times_r = {};
for j = 1:length(conds)
    idx_l = 1;
    idx_r = 1;
    for tr = 1:size(Z.condition, 2)
        if Z.condition(tr) == conds(j) & (Z.TrialError(tr) == 1 | Z.TrialError(tr) == -1) & strcmp(Z.Target_Location{2,tr},'Left')
            event_times_l{conds(j)}(1, idx_l) = PD_Trial(tr);
            event_times_l{conds(j)}(2, idx_l) = PD_Trial(tr) + Z.EyeRT(tr)/1000;
            if(conds(j)~=1 && conds(j)~=3)
                event_times_l{conds(j)}(3, idx_l) = PD_Trial(tr) + Z.ArmRT(tr)/1000;
            end
            idx_l = idx_l + 1;
            %         events_l =
        elseif Z.condition(tr) == conds(j) & (Z.TrialError(tr) == 1 | Z.TrialError(tr) == -1) & strcmp(Z.Target_Location{2,tr},'Right')
            event_times_r{conds(j)}(1, idx_r) = PD_Trial(tr);
            event_times_r{conds(j)}(2, idx_r) = PD_Trial(tr) + Z.EyeRT(tr)/1000;
            if(conds(j)~=1 && conds(j)~=3)
                event_times_r{conds(j)}(3, idx_r) = PD_Trial(tr) + Z.ArmRT(tr)/1000;
            end
            idx_r = idx_r + 1;
        end
    end
end


cd(Folders.KiloFolder)

folderName = 'myFolder';  % Replace with your desired folder name

if ~exist(ses_name_main, 'dir')
    mkdir(ses_name_main);
    disp(['Folder "', ses_name_main, '" created.']);
else
    disp(['Folder "', ses_name_main, '" already exists.']);
end

cd(ses_name_main)
save('event_times_l', 'event_times_l')
save('event_times_r', 'event_times_r')

conditions = unique(Z.condition);

save('conditions', 'conditions')


end