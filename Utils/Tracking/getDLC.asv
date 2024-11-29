function [] = getDLC(Folders, meta)
% For Be240328, there is one condition at the end that I'm not plotting the
% second PD. I used the first 100 pulses for average diff.


% TODO:
% - calibration -- translate to degree
tic
ses_name = meta.sesName;

cd(Folders.save_dir)
% run DLC analyze video first if you have not
if(~exist([ses_name, '_filtered.csv']))
    cd(Folders.TrackingUtils)
    anaconda_path = 'C:\ProgramData\anaconda3';
    setenv('PATH', [anaconda_path ';' anaconda_path '\Scripts;' anaconda_path '\Library\bin;' getenv('PATH')]);
    command = strcat("activate DEEPLABCUT2 && python analyzeVideo.py ", Folders.TrackingUtils, "Belle_Reach_Tracking_DLC_Project/config.yaml", " ", Folders.save_dir, ses_name, '.mp4');
    [status, cmdout] = system(command, "-echo");
end

cd(Folders.save_dir)

if(~exist([ses_name, '_filtered.csv']))
    fileList = dir('*filtered.csv');
    oldFileName = fileList(1).name;  % Assumes only one match; modify if needed
    oldFilePath = fullfile('.', oldFileName);
    % Specify the new name for the file
    newFileName = [ses_name, '_filtered.csv'];  % Change this to your desired name
    newFilePath = fullfile('.', newFileName);
    
    % Rename the file
    movefile(oldFilePath, newFilePath);
end

%% Load PKL file of timestamps
cd(Folders.TrackingUtils)
py.convertPKLtoMAT.convertPKLtoMAT(Folders.save_dir, ses_name);

cd(Folders.save_dir)
load([ses_name, '_timestamps.mat']) % frames_timestamps

% I couldn't convert the returned values from list to double so I save the
% data using scipy to .mat and load it again in the next section

% I think we should restart Matlab each time you change the python code

load(['Z_', ses_name, '_Analog.mat'])
load(['Z_', ses_name, '_Digital.mat'])
load(['Z_', ses_name, '_ML.mat'])

%% Import DLC data

filenamecsv = [ses_name, '_filtered.csv']; 
data_struct = importdata(filenamecsv); 
rawdata = data_struct.data; % import data as matrix in a cell
DLC_filt = rawdata(:, [2, 3, 5, 6, 8, 9, 11, 12, 14, 15])';

ZERO = [416, 156];
DLC_filt(1:2:end, :) = -(DLC_filt(1:2:end, :) - ZERO(1));
DLC_filt(2:2:end, :) = -(DLC_filt(2:2:end, :) - ZERO(2));

% Now check if they have the same size:
if size(DLC_filt, 2) == length(frames_timestamps)
    disp('DLC and timestamp sizes matched!')
else
    disp('DLC and timestamp sizes DO NOT matched!')
end

%% Find Photodiode Values based on frames

if(~exist('ledValue.mat'))
    disp('Getting led values for PD')
    v = VideoReader([ses_name, '.mp4']); % 1 is the camera number
    
    ledLocation = [3, 181]; % Replace x and y with the actual coordinates of PD if you change the framing
    % 3, 320 --> for initial videos
%    v.CurrentTime = 3120;
    ledValue = [];
    frameNum = 1;
    while hasFrame(v)
    
        videoFrame = readFrame(v);
    
        grayFrame = rgb2gray(videoFrame);
    
%         imshow(grayFrame)
%         title(frameNum)
%         pause(0.1)
%         v.CurrentTime
        ledValue = [ledValue, grayFrame(ledLocation(2), ledLocation(1))];
        frameNum = frameNum + 1;
    end
    
    save(['ledValue'], 'ledValue')
else
    load('ledValue')
end
% 
% ledValue(find(frames_timestamps>1900)) = [];
% frames_timestamps(find(frames_timestamps>1900)) = [];

%% Find and plot the starting points

threshold = 10;
data = ledValue;
rising_edges = [];
min_duration = 25; % Minimum number of samples for a valid drop

% Important: In delayed task, we have two or three pulses coming right
% after each other, here, the goal is to just align the two signals, so I
% do not care about every single pulse. For working with the data later,
% you need to use the timepoints from ripple signal and then find the
% corresponding timeopints on the frame signals.

% Assume initial state based on first data point being above or below the threshold
if data(1) > threshold
    lookingForRising = false; % If we start above the threshold, look for a falling edge first
else
    lookingForRising = true; % Otherwise, look for a rising edge
end

for i = 2:length(data)
    if lookingForRising
        if data(i) > threshold && data(i-1) <= threshold % Rising edge detected
            rising_edges = [rising_edges, i];
            lookingForRising = false; % Switch to looking for falling edge
        end
    else
        if data(i) <= threshold && data(i-1) > threshold % Potential falling edge detected
            % Check if this is a brief dip (noise)
            valid_fall = true;
            for j = i+1:min(i+min_duration, length(data))
                if data(j) > threshold % If it rises back quickly, it's noise
                    valid_fall = false;
                    break;
                end
            end
            
            if valid_fall
                lookingForRising = true; % Valid falling edge, continue looking for rising
            end
        end
    end
end

figure

plot(ledValue)
hold on
plot(rising_edges, ledValue(rising_edges), 'o')

%% Read PD data recorded using Ripple

Analog_Data_mV = Z_Analog.data(4, :);
Analog_Data_Time_s = (1:length(Analog_Data_mV))/1000;

% I'm aligning data using the average difference between pulse onsets
time_PD1_ana = Z_Analog.PD_sec(find(~isnan(Z_Analog.PD_sec(:, 1))), 1); % first rising edge in analog data



% % Here's the code you need for sessions that you have PD at the start of each
% % trial
% time_PD1_ana = Z_Digital.data(2, find(Z_Digital.data(1, :) == 9))';
% Analog_Data_mV = zeros(1, length(Analog_Data_Time_s));
% for i = 1:length(time_PD1_ana)
%     [~, index] = min(abs(Analog_Data_Time_s - time_PD1_ana(i)));
%     Analog_Data_mV(index) = 1;
% end



time_PD1_cam = frames_timestamps(rising_edges);

% if(mean(time_PD1_ana) < mean(time_PD1_cam))
%     error('You started neural recording before camera recording!')
% end

time_diff = mean(time_PD1_ana(1:100)' - time_PD1_cam(1:100));
frames_timestamps_corr = frames_timestamps+time_diff;

% sample_diff = dsearchn(Analog_Data_Time_s', time_diff);

% Analog_Data_Time_s_corr = Analog_Data_Time_s(sample_diff:end) - Analog_Data_Time_s(sample_diff);
% Analog_Data_mV_corr = Analog_Data_mV(sample_diff:end);

figure
subplot(2, 1, 1)
plot(Analog_Data_Time_s, normalize(Analog_Data_mV,'range',[0,1]))
hold on
plot(frames_timestamps, normalize(ledValue,'range',[0,1]))
title('Un-aligned PDs')

subplot(2, 1, 2)
plot(Analog_Data_Time_s, normalize(Analog_Data_mV,'range',[0,1]))
hold on
plot(frames_timestamps_corr, normalize(ledValue,'range',[0,1]))
title('Aligned PDs')

% There still seems to be a slight difference between the pulses in the 
% Camera and PD signals, possibly due to the camera timestamps
% (which currently doesn't make sense to me).

%% show video for a given trial aligned to PD onset 

% % Now give me trial number, I will plot the hand data for you
% tr = 393; % she reached to the wrong location on tr = 12
% 
% if(Z_ML.condition == 5 | Z_ML.condition == 6)
%     pre_dur_tr = 500; 
% else
%     pre_dur_tr = 200;
% end
% 
% pre_dur_tr = 500; 
% 
% pre_dur = 500; %ms
% post_dur = 500; %ms
% if(~isnan(Z_Analog.PD_sec(tr, 1)))
% 
%     cam_sample_st = dsearchn(frames_timestamps_corr', Z_Analog.PD_sec(tr, 1) - pre_dur/1000);
%     cam_sample_en = dsearchn(frames_timestamps_corr', Z_Analog.PD_sec(tr, 1) + post_dur/1000);
%     
%     v = VideoReader([ses_name, '.mp4']); % 1 is the camera number
%     
%     for frameNumber = cam_sample_st:cam_sample_en
%     
%         % Calculate the frame number based on the frame rate and the specified time
%         time = frameNumber/v.FrameRate;
%         % Seek to the frame directly if using MATLAB R2021a or later
%         v.CurrentTime = time;
%     
%         subplot(5, 2, [1, 3, 5, 7])
% 
%         if hasFrame(v)
%             videoFrame = readFrame(v);
%             img_double = im2double(videoFrame);
%             brighter_img = img_double + 0.2;
%             brighter_img(brighter_img > 1) = 1;
%             brighter_img = imrotate(brighter_img, 180);
%             [N, M, ~] = size(img_double);
% 
%             subplot(5, 2, [1, 3, 5, 7])
%             imshow(brighter_img); % Display the frame
%     
%             hold on
%             x_new = M - (-DLC_filt(1, frameNumber)+ZERO(1)) + 1;
%             y_new = N - (-DLC_filt(2, frameNumber)+ZERO(2)) + 1;
%             plot(x_new, y_new, 'bo')
% 
%             x_new = M - (-DLC_filt(3, frameNumber)+ZERO(1)) + 1;
%             y_new = N - (-DLC_filt(4, frameNumber)+ZERO(2)) + 1;
%             plot(x_new, y_new, 'go')
% 
%             x_new = M - (-DLC_filt(5, frameNumber)+ZERO(1)) + 1;
%             y_new = N - (-DLC_filt(6, frameNumber)+ZERO(2)) + 1;
%             plot(x_new, y_new, 'ko')
% 
%             x_new = M - (-DLC_filt(7, frameNumber)+ZERO(1)) + 1;
%             y_new = N - (-DLC_filt(8, frameNumber)+ZERO(2)) + 1;
%             plot(x_new, y_new, 'mo')
%             
%             subplot(5, 2, 2)
%             plot(1000*(0:1/v.FrameRate:((frameNumber - cam_sample_st)/v.FrameRate)) - pre_dur, DLC_filt(1, cam_sample_st:frameNumber), 'b')
%             title('Dig1')
%             ylim([min(DLC_filt(1, cam_sample_st:cam_sample_en))*1, max(DLC_filt(1, cam_sample_st:cam_sample_en))*1])
%             xline(0, '--')
%             xlim([-pre_dur, post_dur])
%             xline(Z_ML.ArmRT(tr), 'r--')
% 
%             subplot(5, 2, 4)
%             plot(1000*(0:1/v.FrameRate:((frameNumber - cam_sample_st)/v.FrameRate)) - pre_dur, DLC_filt(3, cam_sample_st:frameNumber), 'g')
%             title('Dig2')
%             ylim([min(DLC_filt(3, cam_sample_st:cam_sample_en))*1, max(DLC_filt(3, cam_sample_st:cam_sample_en))*1])
%             xline(0, '--')
%             xlim([-pre_dur, post_dur])
%             xline(Z_ML.ArmRT(tr), 'r--')
% 
%             subplot(5, 2, 6)
%             plot(1000*(0:1/v.FrameRate:((frameNumber - cam_sample_st)/v.FrameRate)) - pre_dur, DLC_filt(5, cam_sample_st:frameNumber), 'k')
%             title('Dig3')
%             ylim([min(DLC_filt(5, cam_sample_st:cam_sample_en))*1, max(DLC_filt(5, cam_sample_st:cam_sample_en))*1])
%             xline(0, '--')
%             xlim([-pre_dur, post_dur])
%             xline(Z_ML.ArmRT(tr), 'r--')
% 
%             subplot(5, 2, 8)
%             plot(1000*(0:1/v.FrameRate:((frameNumber - cam_sample_st)/v.FrameRate)) - pre_dur, DLC_filt(7, cam_sample_st:frameNumber), 'm')
%             title('Dig4')
%             xlabel('time (s)')
%             ylim([min(DLC_filt(7, cam_sample_st:cam_sample_en))*1, max(DLC_filt(7, cam_sample_st:cam_sample_en))*1])
%             xline(0, '--')
%             xlim([-pre_dur, post_dur])
%             xline(Z_ML.ArmRT(tr), 'r--')
% 
%             time_vector_fifth = -pre_dur_tr:(length(Z_ML.TouchX_Diode{20}) - pre_dur_tr);
%             time_frame = round(1000* (time - cam_sample_st/v.FrameRate) - pre_dur);
%             if time_frame>= -pre_dur_tr
%                 idx_in_range = find(time_vector_fifth <= time_frame, 1, 'last'); % Find the closest index up to the current time
%                 if ~isempty(idx_in_range)
%                     subplot(5, 2, 10)
%                     plot(time_vector_fifth(1:idx_in_range), Z_ML.TouchX_Diode{tr}(1:idx_in_range))
%                     title('ML Touchscreen')
%                     xlabel('time (ms)')
%                     ylim([min(Z_ML.TouchX_Diode{tr})*1, max(Z_ML.TouchX_Diode{tr})*1])
%                 end
%             end
%             xline(0, '--')
%             xlim([-pre_dur, post_dur])
%             xline(Z_ML.ArmRT(tr), 'r--')
% 
%             sgtitle(['Trial number: ', num2str(tr), ',   ', num2str(round(1000* (time - cam_sample_st/v.FrameRate)))])
% 
%         end
%     end
% 
% end

%% Save data

Z_TR.data = DLC_filt;
Z_TR.frames_timestamps = frames_timestamps_corr;
Z_TR.PD_px_Value = ledValue;

save(['Z_', ses_name, '_TR'], 'Z_TR')
disp('Z_TR is saved!')

toc
end