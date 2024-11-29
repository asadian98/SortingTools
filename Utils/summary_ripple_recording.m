function [] = summary_ripple_recording(animal_name, Researcher, ses_name, tail)

% This function reads the ripple files and gives you a summary of your recording.
% There is no label for raw and LFP channels showing which electrode is
% that. So you need to know the channel numbers yourself.

% animal_name = 'Belle';
% Researcher = 'test_KS';
% ses_name = 'Be240404';
% RP_tail_list = {'_001'};
% tail = '_001';

Folders.O_NHP_Folder = ['C:\Users\CorneilLab\Desktop\', Researcher, '\'];
Folders.NHP_Folder = animal_name;
Folders.Ripple_Raw_Data_Folder = [Folders.O_NHP_Folder,'Ripple_Raw_Data\'];
Folders.ML_Behav_Folder = [Folders.O_NHP_Folder,'\ML_Behavioural_Files\', animal_name];
Folders.KiloFolder = [Folders.O_NHP_Folder, Folders.NHP_Folder, '\'];
Folders.save_dir = [Folders.O_NHP_Folder, Folders.NHP_Folder, '\', ses_name];

completeFilePath_Digital = [Folders.Ripple_Raw_Data_Folder,animal_name, '\', [ses_name, tail],'.nev'];

% Open the file and extract some basic information
% disp('opening Digital datafile to acquire statecodes')
[ns_status, hFile] = ns_OpenFile(completeFilePath_Digital,'single');  % opens only .nev file
if strcmp(ns_status,'ns_OK')
    disp('NEV opened properly')
else
    disp('Did not open NEV properly')
end

% Assuming hFile.Entity(:).Label contains the labels

labels = {hFile.Entity(:).Label};  % Extract labels into a cell array
entity_types = {hFile.Entity(:).EntityType};  % Extract entity types

num_labels = length(labels);       % Get the number of labels

% Initialize a map to store electrode contact counts and indices
electrode_data = containers.Map();
hasEvents = 0;

for i = 1:num_labels

    % Check if the EntityType is 'Event' (statecodes)
    if strcmp(entity_types{i}, 'Event')
        hasEvents = 1;  % Skip entities of type 'Event'
    end

    label = labels{i};

    % Skip empty or invalid labels (such as {0x0 char})
    if isempty(label) || (iscell(label) && isempty(label{1}))
        continue;
    end

    % Extract the electrode part and the contact number, ignore 'spk'
    parts = strsplit(label, '-');
    if length(parts) >= 3
        % Extract the electrode part and contact number
        electrode = strcat(parts{1}, '-', parts{2});   % Electrode part (e.g., 'A-1')
        contact_number = strsplit(parts{3}, ' ');      % Split at space to remove 'spk'
        contact_number = contact_number{1};            % Get the contact number without 'spk'
    else
        continue;  % If the label format is unexpected, skip it
    end

    % Update the count and store the index and contact number for this electrode
    if isKey(electrode_data, electrode)
        % Increment the count and add the index and contact number
        current_data = electrode_data(electrode);
        current_data.count = current_data.count + 1;
        current_data.indices = [current_data.indices, i];           % Add current index
        current_data.contact_numbers = [current_data.contact_numbers, {contact_number}];  % Add contact number as cell array
        electrode_data(electrode) = current_data;  % Update the map entry
    else
        % Create a new entry with count 1 and store the current index and contact number
        electrode_data(electrode) = struct('count', 1, 'indices', i, 'contact_numbers', {contact_number});
    end
end

% Display the summary with counts, indices, and contact numbers
keys = electrode_data.keys();
for i = 1:length(keys)
    electrode = keys{i};
    data = electrode_data(electrode);
    fprintf('Spikes: electrode %s has %d contacts with contact numbers: %s\n',...
        electrode, data.count, strjoin(data.contact_numbers, ', '));
end

if(hasEvents)
    disp('The recording includes statecodes.')
else
    disp('The recording does NOT include statecodes.')
end

%% LFP and Analog data

completeFilePath_Digital = [Folders.Ripple_Raw_Data_Folder,animal_name, '\', [ses_name, tail],'.ns2'];

% Open the file and extract some basic information
[ns_status, hFile] = ns_OpenFile(completeFilePath_Digital,'single');  % opens only .nev file
if strcmp(ns_status,'ns_OK')
    disp('NS2 opened properly')
else
    disp('Did not open NS2 properly')
end

% Assuming hFile.Entity(:).Label contains the labels

labels = {hFile.Entity(:).Label};  % Extract labels into a cell array
entity_types = {hFile.Entity(:).EntityType};  % Extract entity types

num_labels = length(labels);       % Get the number of labels

% Initialize counters for 'lfp' and 'analog'
lfp_count = 0;
analog_count = 0;

% Loop through the labels and count how many are 'lfp' and 'analog'
for i = 1:num_labels
    label = labels{i};

    % Check if label starts with 'lfp'
    if startsWith(label, 'lfp')
        lfp_count = lfp_count + 1;
        % Check if label starts with 'analog'
    elseif startsWith(label, 'analog')
        analog_count = analog_count + 1;
    end
end

% Display the results
fprintf('Number of LFP contacts: %d\n', lfp_count);
fprintf('Number of Analog contacts: %d\n', analog_count);

%% Raw Files

completeFilePath_Digital = [Folders.Ripple_Raw_Data_Folder,animal_name, '\', [ses_name, tail],'.ns5'];

% Open the file and extract some basic information
[ns_status, hFile] = ns_OpenFile(completeFilePath_Digital,'single');  % opens only .nev file
if strcmp(ns_status,'ns_OK')
    disp('NS5 opened properly')
else
    disp('Did not open NS5 properly')
end

% Assuming hFile.Entity(:).Label contains the labels

labels = {hFile.Entity(:).Label};  % Extract labels into a cell array
entity_types = {hFile.Entity(:).EntityType};  % Extract entity types

num_labels = length(labels);       % Get the number of labels

% Initialize counters for 'raw'
raw_count = 0;

% Loop through the labels and count how many are 'raw'
for i = 1:num_labels
    label = labels{i};

    % Check if label starts with 'lfp'
    if startsWith(label, 'raw')
        raw_count = raw_count + 1;
    end
end

% Display the results
fprintf('Number of Raw contacts: %d\n', raw_count);

end