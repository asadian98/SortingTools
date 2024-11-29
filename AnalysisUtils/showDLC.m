%% show video for a given trial aligned to PD onset 

ses_name = 'Be240423';
cd(['C:\Users\CorneilLab\Desktop\AA\Belle\', ses_name])
load(['Z_', ses_name, '_ML.mat'])
load(['Z_', ses_name, '_TR.mat'])
load(['Z_', ses_name, '_Analog.mat'])
load(['Z_', ses_name, '_Digital.mat'])
ledValue = Z_TR.PD_px_Value;

% Get pulse onsets first and compare it with Z_Analog.PD_sec

threshold = 10;
data = ledValue;
rising_edges = [];
min_duration = 25; % Minimum number of samples for a valid drop

% Assume initial state based on first data point being above or below the threshold
if data(1) > threshold
    lookingForRising = false; % If we start above the threshold, look for a falling edge first
else
    lookingForRising = true; % Otherwise, look for a rising edge
end

for i = 2:length(data)
    if lookingForRising
        if data(i) > threshold && data(i-1) <= threshold % Rising edge detected
            rising_edges = [rising_edges, i-1];
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

if(length(rising_edges) == length(find(~isnan(Z_Analog.PD_sec(:, 1)))))
    disp('You are good to go')
else
    disp('Sizes do not match!')
end

DLC_filt = Z_TR.data;
ZERO = [416, 156];

%%

% diff_RT = [];
% tr_arr = [];
% Z_ML.ArmRT_DLC = zeros(1, length(Z_ML.condition)) - 999;
% 
% for tr = 1:length(Z_ML.condition)
% 
%     if(Z_ML.condition(tr) == 5 | Z_ML.condition(tr) == 6)
%         pre_dur_tr = 500; 
%     else
%         pre_dur_tr = 200;
%     end
%     pre_dur = 0; %ms
%     post_dur = 500; %ms
%     vel_thresh = 0.2; % in percentage
% 
%     if(~isnan(Z_Analog.PD_sec(tr, 1)) & Z_ML.ArmRT(tr) > 0)
%     
%         PD_onset_ripple_in_DLC = dsearchn(Z_TR.frames_timestamps', Z_Analog.PD_sec(tr, 1));
%         cam_sample_st = dsearchn(Z_TR.frames_timestamps', Z_TR.frames_timestamps(rising_edges(dsearchn(rising_edges', PD_onset_ripple_in_DLC))) - pre_dur/1000);
%         cam_sample_en = dsearchn(Z_TR.frames_timestamps', Z_TR.frames_timestamps(rising_edges(dsearchn(rising_edges', PD_onset_ripple_in_DLC))) + post_dur/1000);
%     
%         % Find DLC movement onset
%         onset_time_min = 10000;
%         for dig_idx = 1:4
%             position = DLC_filt(dig_idx, cam_sample_st:cam_sample_en);
%             % Calculate velocity if not already provided
%             velocity = diff(position); % Simple difference method
%             velocity = [0, velocity]; % Append 0 to match the length with 'position'
%             velocity_smooth = smoothdata(velocity, 'movmean', 5);
%             % Calculate the maximum speed and set the threshold to 10% of it
%             max_speed = max(abs(velocity_smooth)); % Maximum speed (absolute value)
%             threshold = vel_thresh * max_speed; % 10% of maximum speed
%             % Find the movement onset index
%             onset_index = find(abs(velocity_smooth) > threshold, 1, 'first');
%             DLC_time_Array = 1000*(0:1/v.FrameRate:((cam_sample_en - cam_sample_st)/v.FrameRate)) - pre_dur;
%             onset_time = DLC_time_Array(onset_index);
%             if(onset_time < onset_time_min)
%                 onset_time_min = onset_time;
%             end
%         end
% 
%         tr_arr = [tr_arr, tr];
%         diff_RT = [diff_RT, Z_ML.ArmRT(tr) - onset_time_min];
%         Z_ML.ArmRT_DLC(tr) = onset_time_min;
% 
%     end
% end
% 
% disp(['Average diff  =  ', num2str(mean(diff_RT))])
% save(['Z_', ses_name, '_ML.mat'], 'Z_ML')

%%

% Now give me trial number, I will plot the hand data for you
tr = 492; % she reached to the wrong location on tr = 12

% time_PD1_ana = Z_Analog.PD_sec(tr, 1);

if(Z_ML.condition(tr) == 5 | Z_ML.condition(tr) == 6)
    pre_dur_tr = 500; 
else
    pre_dur_tr = 200;
end

% pre_dur_tr = 500; 

pre_dur = 500; %ms
post_dur = 500; %ms
vel_thresh = 0.2; % in percentage

% Create a new figure
fig = figure;
set(fig, 'Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.9]); % Adjusts to 90% of screen size with small margins

if(~isnan(Z_Analog.PD_sec(tr, 1)))

    PD_onset_ripple_in_DLC = dsearchn(Z_TR.frames_timestamps', Z_Analog.PD_sec(tr, 1));
    cam_sample_st = dsearchn(Z_TR.frames_timestamps', Z_TR.frames_timestamps(rising_edges(dsearchn(rising_edges', PD_onset_ripple_in_DLC))) - pre_dur/1000);
    cam_sample_en = dsearchn(Z_TR.frames_timestamps', Z_TR.frames_timestamps(rising_edges(dsearchn(rising_edges', PD_onset_ripple_in_DLC))) + post_dur/1000);
    
    onset_time_min = Z_ML.ArmRT_DLC(tr);
    v = VideoReader([ses_name, '.mp4']); % 1 is the camera number
    
    % Set up the video writer object
    video = VideoWriter(['DLC_tr', num2str(tr), '.avi']); % Specify the file name and format
    video.FrameRate = v.FrameRate; % Set the frame rate
    open(video); % Open the video file for writing

    for frameNumber = cam_sample_st:cam_sample_en
        DLC_time_Array = 1000*(0:1/v.FrameRate:((frameNumber - cam_sample_st)/v.FrameRate)) - pre_dur;

        % Calculate the frame number based on the frame rate and the specified time
        time = frameNumber/v.FrameRate;
        % Seek to the frame directly if using MATLAB R2021a or later
        v.CurrentTime = time;
    
        subplot(5, 2, [1, 3, 5, 7])

        if hasFrame(v)
            videoFrame = readFrame(v);
            img_double = im2double(videoFrame);
            brighter_img = img_double + 0.2;
            brighter_img(brighter_img > 1) = 1;
            brighter_img = imrotate(brighter_img, 180);
            [N, M, ~] = size(img_double);

            subplot(5, 2, [1, 3, 5, 7])
            imshow(brighter_img); % Display the frame
    
            hold on
            x_new = M - (-DLC_filt(1, frameNumber)+ZERO(1)) + 1;
            y_new = N - (-DLC_filt(2, frameNumber)+ZERO(2)) + 1;
            plot(x_new, y_new, 'bo')

            x_new = M - (-DLC_filt(3, frameNumber)+ZERO(1)) + 1;
            y_new = N - (-DLC_filt(4, frameNumber)+ZERO(2)) + 1;
            plot(x_new, y_new, 'go')

            x_new = M - (-DLC_filt(5, frameNumber)+ZERO(1)) + 1;
            y_new = N - (-DLC_filt(6, frameNumber)+ZERO(2)) + 1;
            plot(x_new, y_new, 'ko')

            x_new = M - (-DLC_filt(7, frameNumber)+ZERO(1)) + 1;
            y_new = N - (-DLC_filt(8, frameNumber)+ZERO(2)) + 1;
            plot(x_new, y_new, 'mo')
            
            subplot(5, 2, 2)
            plot(DLC_time_Array, DLC_filt(1, cam_sample_st:frameNumber), 'b')
            title('Dig1')
            ylim([min(DLC_filt(1, cam_sample_st:cam_sample_en))*1, max(DLC_filt(1, cam_sample_st:cam_sample_en))*1])
            xline(0, '--')
            xlim([-pre_dur, post_dur])
            xline(Z_ML.ArmRT(tr), 'r--')
            xline(onset_time_min, 'b--')

% If you wanna plot the PD on DLC
%             subplot(5, 2, 2)
%             plot(DLC_time_Array, ledValue(1, cam_sample_st:frameNumber), 'b')
%             title('Dig1')
% %             ylim([min(ledValue(1, cam_sample_st:cam_sample_en))*1, max(ledValue(1, cam_sample_st:cam_sample_en))*1])
%             xline(0, '--')
%             xlim([-pre_dur, post_dur])
%             xline(Z_ML.ArmRT(tr), 'r--')
%             xline(onset_time_min, 'b--')

            subplot(5, 2, 4)
            plot(DLC_time_Array, DLC_filt(3, cam_sample_st:frameNumber), 'g')
            title('Dig2')
            ylim([min(DLC_filt(3, cam_sample_st:cam_sample_en))*1, max(DLC_filt(3, cam_sample_st:cam_sample_en))*1])
            xline(0, '--')
            xlim([-pre_dur, post_dur])
            xline(Z_ML.ArmRT(tr), 'r--')
            xline(onset_time_min, 'b--')

            subplot(5, 2, 6)
            plot(DLC_time_Array, DLC_filt(5, cam_sample_st:frameNumber), 'k')
            title('Dig3')
            ylim([min(DLC_filt(5, cam_sample_st:cam_sample_en))*1, max(DLC_filt(5, cam_sample_st:cam_sample_en))*1])
            xline(0, '--')
            xlim([-pre_dur, post_dur])
            xline(Z_ML.ArmRT(tr), 'r--')
            xline(onset_time_min, 'b--')

            subplot(5, 2, 8)
            plot(DLC_time_Array, DLC_filt(7, cam_sample_st:frameNumber), 'm')
            title('Dig4')
            xlabel('time (s)')
            ylim([min(DLC_filt(7, cam_sample_st:cam_sample_en))*1, max(DLC_filt(7, cam_sample_st:cam_sample_en))*1])
            xline(0, '--')
            xlim([-pre_dur, post_dur])
            xline(Z_ML.ArmRT(tr), 'r--')
            xline(onset_time_min, 'b--')

            time_vector_fifth = -pre_dur_tr:(length(Z_ML.TouchX_Diode{tr}) - pre_dur_tr);
            time_frame = round(1000* (time - cam_sample_st/v.FrameRate) - pre_dur);
            if time_frame>= -pre_dur_tr
                idx_in_range = find(time_vector_fifth <= time_frame, 1, 'last'); % Find the closest index up to the current time
                if ~isempty(idx_in_range)
                    subplot(5, 2, 10)
                    plot(time_vector_fifth(1:idx_in_range), Z_ML.TouchX_Diode{tr}(1:idx_in_range))
                    title('ML Touchscreen')
                    xlabel('time (ms)')
                    ylim([min(Z_ML.TouchX_Diode{tr})*1, max(Z_ML.TouchX_Diode{tr})*1])
                    ylim([-20, 20])
                end
            end
            xline(0, '--')
            xlim([-pre_dur, post_dur])
            xline(Z_ML.ArmRT(tr), 'r--')
            xline(onset_time_min, 'b--')

            sgtitle(['Trial number: ', num2str(tr), ',   ', num2str(round(1000* (time - cam_sample_st/v.FrameRate)))])
            
            % Capture the current figure as a frame
            frame = getframe(gcf);
            
            % Write the frame to the video
            writeVideo(video, frame);
        end
    end

end

% Close the video file
close(video);

disp('Animation saved successfully!');