function analyze_NHP_plotAll(Z_ML, Z_Analog, Z_Spikes, ch_want, cond_want, alignment)

%   analyze_NHP_plotAll.m Plots all units for left and right movements
%
%   Inputs:
%       Z Files,
%       ch_want: the channel number you want to look at
%       cond_want: the condition you want
%       alignment: 'Target', 'Saccade', 'Arm'
%
%   Outputs:
%       no output, it will generate a figure and plot your data 
%       
%   TODO: alignment for the delayed task, make ROC working
%   

lowRT_cutoff = 70;
highRT_cutoff = 300;

f = figure;
panel1 = uipanel('Parent',f);
panel2 = uipanel('Parent',panel1);
set(panel1,'Position',[0 0 0.95 1]);
set(panel2,'Position',[0 -3 1 4]);
% h = image;
set(gca,'Parent',panel2);
s = uicontrol('Style','Slider','Parent',f,...
    'Units','normalized','Position',[0.95 0 0.05 1],...
    'Value',1,'Callback',{@slider_callback1,panel2});
t = tiledlayout(panel2, length(ch_want), max([Z_Spikes.data.nUnits]));
t.TileSpacing = "tight";
t.Padding = "tight"; % tight, none, compact

st = 500;
en = 500;
time_re_target = -st:en;

LEN = en+st+1;
flagg = 0;
r = 0;
l = [];

for ch = ch_want

    channelFlag = 0;

    for unitID = 1:Z_Spikes.data(ch).nUnits

        trs_right = find(Z_ML.condition == cond_want & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Right'));
        trs_left = find(Z_ML.condition == cond_want & Z_ML.TrialError == 1 & strcmp(Z_ML.Target_Location(2,:),'Left'));
        spk = Z_Spikes.data(ch).neural_timeStamps(1, find(Z_Spikes.data(ch).neural_timeStamps(2, :) == Z_Spikes.data(ch).Units(unitID)));

        [SPIKES_RIGHT, SPDEN_RIGHT] = get_spikes_trial(Z_Analog, spk, trs_right);
        [SPIKES_LEFT, SPDEN_LEFT]  = get_spikes_trial(Z_Analog, spk, trs_left);

        ArmON_LEFT = Z_ML.ArmRT(trs_left);
        ArmON_RIGHT = Z_ML.ArmRT(trs_right);
        EYE_RT_LEFT = Z_ML.EyeRT(trs_left);
        EYE_RT_RIGHT = Z_ML.EyeRT(trs_right);


        % If there are spikes between 40 and 100 ms on 10% of all trials to the
        % left and right, then proceed with plotting this figure
        [rows_RIGHT,columns_RIGHT]=find(SPIKES_RIGHT >= 40 & SPIKES_RIGHT <=300);
        rows_RIGHT = unique(rows_RIGHT);   % Row identifies with at least one spike between 40 and 100 ms
        [rows_LEFT,columns_LEFT]=find(SPIKES_LEFT >= 40 & SPIKES_LEFT <=300);
        rows_LEFT = unique(rows_LEFT);   % Row identifies with at least one spike between 40 and 100 ms


        if length(rows_LEFT) + length(rows_RIGHT) >= round((length(EYE_RT_LEFT) + length(EYE_RT_RIGHT))/100) && size(SPDEN_RIGHT, 1) > 1 && (size(SPDEN_LEFT, 1) > 1 || cond_want == 7)

            % ALIGN DATA ON RT
            pre_RT = -st;
            post_RT = en;
            time_re_sac = pre_RT:post_RT;
            SPIKES_RIGHT_RT = [];
            SPDEN_RIGHT_RT = zeros(length(EYE_RT_RIGHT),length(time_re_sac));
            SPIKES_LEFT_RT = [];
            SPDEN_LEFT_RT = zeros(length(EYE_RT_LEFT),length(time_re_sac));
            for i = 1:size(SPDEN_RIGHT, 1)
                SPIKES_RIGHT_RT(i,:) = SPIKES_RIGHT(i,:)-EYE_RT_RIGHT(i); % Shift raster timepoints by RT

                RTpoint = find(time_re_target == round(EYE_RT_RIGHT(i)));
                RT_pre_cut = RTpoint+pre_RT; RT_post_cut = RTpoint+post_RT;
                if RT_post_cut > length(time_re_target); RT_post_cut = length(time_re_target); end
                snippet = SPDEN_RIGHT(i,RT_pre_cut:RT_post_cut);  % This is the snippet around RT
                SPDEN_RIGHT_RT(i,1:length(snippet)) = snippet;
            end

            if(cond_want ~= 7)
                for i = 1:size(SPDEN_LEFT, 1)
                    SPIKES_LEFT_RT(i,:) = SPIKES_LEFT(i,:)-EYE_RT_LEFT(i); % Shift raster timepoints by RT

                    RTpoint = find(time_re_target == round(EYE_RT_LEFT(i)));
                    RT_pre_cut = RTpoint+pre_RT; RT_post_cut = RTpoint+post_RT;
                    if RT_post_cut > length(time_re_target); RT_post_cut = length(time_re_target); end
                    snippet = SPDEN_LEFT(i,RT_pre_cut:RT_post_cut);  % This is the snippet around RT
                    SPDEN_LEFT_RT(i,1:length(snippet)) = snippet;
                end
            end

            % ALIGN DATA ON ARM RT
            pre_RT = -st;
            post_RT = en;

            time_re_sac = pre_RT:post_RT;
            SPIKES_RIGHT_ARM_RT = [];
            SPDEN_RIGHT_ARM_RT = zeros(length(ArmON_RIGHT),length(time_re_sac));
            SPIKES_LEFT_ARM_RT = [];
            SPDEN_LEFT_ARM_RT = zeros(length(ArmON_LEFT),length(time_re_sac));
            for i = 1:size(SPDEN_RIGHT, 1)
                SPIKES_RIGHT_ARM_RT(i,:) = SPIKES_RIGHT(i,:)-ArmON_RIGHT(i); % Shift raster timepoints by RT
                RTpoint = find(time_re_target == round(ArmON_RIGHT(i)));
                RT_pre_cut = RTpoint+pre_RT; RT_post_cut = RTpoint+post_RT;
                if RT_post_cut > length(time_re_target); RT_post_cut = length(time_re_target); end
                snippet = SPDEN_RIGHT(i,RT_pre_cut:RT_post_cut);  % This is the snippet around RT
                SPDEN_RIGHT_ARM_RT(i,1:length(snippet)) = snippet;
            end

            for i = 1:size(SPDEN_LEFT, 1)
                %             size(SPIKES_LEFT(i,:))
                %             size(ArmON_LEFT(i))
                SPIKES_LEFT_ARM_RT(i,:) = SPIKES_LEFT(i,:)-ArmON_LEFT(i); % Shift raster timepoints by RT

                RTpoint = find(time_re_target == round(ArmON_LEFT(i)));
                RT_pre_cut = RTpoint+pre_RT; RT_post_cut = RTpoint+post_RT;
                if RT_post_cut > length(time_re_target); RT_post_cut = length(time_re_target); end
                snippet = SPDEN_LEFT(i,RT_pre_cut:RT_post_cut);  % This is the snippet around RT
                SPDEN_LEFT_ARM_RT(i,1:length(snippet)) = snippet;
            end

            % ALIGN DATA ON Reach GO
            if(cond_want == 8)
                pre_RT = -st;
                post_RT = en;
                time_re_sac = pre_RT:post_RT;
                SPIKES_RIGHT_RGO_RT = [];
                SPDEN_RIGHT_RGO_RT = zeros(size(SPDEN_RIGHT, 1),length(time_re_sac));
                SPIKES_LEFT_RGO_RT = [];
                SPDEN_LEFT_RGO_RT = zeros(size(SPDEN_LEFT, 1),length(time_re_sac));
                for i = 1:size(SPDEN_RIGHT, 1)
                    SPIKES_RIGHT_RGO_RT(i,:) = SPIKES_RIGHT(i,:)-ReachPhotoOffset_RIGHT(i); % Shift raster timepoints by RT
                    RTpoint = find(time_re_target == round(ReachPhotoOffset_RIGHT(i)));
                    RT_pre_cut = RTpoint+pre_RT; RT_post_cut = RTpoint+post_RT;
                    if RT_post_cut > length(time_re_target); RT_post_cut = length(time_re_target); end
                    snippet = SPDEN_RIGHT(i,RT_pre_cut:RT_post_cut);  % This is the snippet around RT
                    SPDEN_RIGHT_RGO_RT(i,1:length(snippet)) = snippet;
                end

                for i = 1:size(SPDEN_LEFT, 1)
                    %             size(SPIKES_LEFT(i,:))
                    %             size(ArmON_LEFT(i))
                    SPIKES_LEFT_RGO_RT(i,:) = SPIKES_LEFT(i,:)-ReachPhotoOffset_LEFT(i); % Shift raster timepoints by RT
                    RTpoint = find(time_re_target == round(ReachPhotoOffset_LEFT(i)));
                    RT_pre_cut = RTpoint+pre_RT; RT_post_cut = RTpoint+post_RT;
                    if RT_post_cut > length(time_re_target); RT_post_cut = length(time_re_target); end
                    snippet = SPDEN_LEFT(i,RT_pre_cut:RT_post_cut);  % This is the snippet around RT
                    SPDEN_LEFT_RGO_RT(i,1:length(snippet)) = snippet;
                end
            end

            % Now sort data based on EYE RT
            [EYE_RT_LEFT_sorted,I_EYE_RT_LEFT_sorted]=sort(EYE_RT_LEFT);
            [EYE_RT_RIGHT_sorted,I_EYE_RT_RIGHT_sorted]=sort(EYE_RT_RIGHT);

            % Now sort data based on ARM RT
            [ARM_RT_LEFT_sorted,I_ARM_RT_LEFT_sorted]=sort(ArmON_LEFT);
            [ARM_RT_RIGHT_sorted,I_ARM_RT_RIGHT_sorted]=sort(ArmON_RIGHT);

            SPIKES_LEFT_sorted = SPIKES_LEFT(I_EYE_RT_LEFT_sorted,:);
            SPIKES_RIGHT_sorted = SPIKES_RIGHT(I_EYE_RT_RIGHT_sorted,:);
            SPIKES_LEFT_RT_sorted = SPIKES_LEFT_RT(I_EYE_RT_LEFT_sorted,:);
            SPIKES_RIGHT_RT_sorted = SPIKES_RIGHT_RT(I_EYE_RT_RIGHT_sorted,:);
            SPIKES_LEFT_ARM_RT_sorted = SPIKES_LEFT_ARM_RT(I_ARM_RT_LEFT_sorted,:);
            SPIKES_RIGHT_ARM_RT_sorted = SPIKES_RIGHT_ARM_RT(I_ARM_RT_RIGHT_sorted,:);

            nexttile(t, (find(ch_want == ch)-1)*max([Z_Spikes.data.nUnits])+unitID)

            if(channelFlag == 0)
                ylabel(['Ch #', num2str(find(ch_want == ch))])
                channelFlag = 1;
            end

            if(strcmp(alignment, 'Target'))
                R_data = SPDEN_RIGHT;
                L_data = SPDEN_LEFT;
            elseif(strcmp(alignment, 'Saccade'))
                R_data = SPDEN_RIGHT_RT;
                L_data = SPDEN_LEFT_RT;
            elseif(strcmp(alignment, 'Arm'))
                R_data = SPDEN_RIGHT_ARM_RT;
                L_data = SPDEN_LEFT_ARM_RT;
            elseif(strcmp(alignment, 'ReachGO'))
                R_data = SPDEN_RIGHT_RGO_RT;
                L_data = SPDEN_LEFT_RGO_RT;
            end

            plot_patchplot(-st:en,R_data,'r', '');
            if(cond_want ~= 7)
                plot_patchplot(-st:en,L_data,'g','');
            end
            xlim([-st, en])
            title(['Unit#', num2str(unitID)])
            xline(0, '--')

            % RUN ROC ANALYSIS ON SPDEN FUNCTIONS.
            min_ROC_time = -100; max_ROC_time = 500;
            ROC_time = min_ROC_time: max_ROC_time;
            ROC = [];

            %                     tic
            if(cond_want ~= 7)
                for i = 1:length(ROC_time)
                    colmatch = find(-st:en == ROC_time(i));
                    % col1: time col2: AUC Value
                    ROC(end+1,:) = [ROC_time(i), calcROC_BC(SPDEN_RIGHT(:,colmatch), SPDEN_LEFT(:,colmatch))];
                end
                %                     toc

                %                     subplot(gridV,3,4); hold on; title('Time series ROC on spike density functions')
                %                     plot(ROC(:,1),ROC(:,2),'b-')
                %                     axis([min_ROC_time max_ROC_time 0 1])
                %                     plot([0 0],[0 1],'k--')
                %                     plot([min_ROC_time max_ROC_time],[.5 .5],'k--')
                %                     plot([min_ROC_time max_ROC_time],[.6 .6],'k-.')
                %                     plot([min_ROC_time max_ROC_time],[.4 .4],'k-.')
                onset_latency = find_onset_latency_ROC(ROC_time,ROC(:,2),[30 100],[.6 8 10],0);
            else
                onset_latency = 0;
            end
            if isnan(onset_latency)
                onset_latency = 0;
                %                         text(min_ROC_time + 50,0.1,strcat('ROC discrmination time =',num2str(onset_latency)));
                %                         xline(onset_latency, 'k-.')
            end

            baseline_ROC_start_index = find(ROC_time == -100);
            baseline_ROC_end_index = find(ROC_time == 30);
            if onset_latency ~= 0 & sum(ROC(baseline_ROC_start_index:baseline_ROC_end_index,2)) > 0
                % Then detrend the ROC from -100 to +30, replot, and
                % recalculate threshold
                R = linfit_jo(-100:30,ROC(baseline_ROC_start_index:baseline_ROC_end_index,2)');
                lineartrend = [];
                for x = min_ROC_time:max_ROC_time
                    lineartrend(end+1) = R(1)*x+R(3)-0.5;   % to center it around 0.5
                end
                ROC_detrend = ROC(:,2)-lineartrend';

                %                         subplot(gridV,3,5); hold on; title('Detrended time series ROC on spike density functions')
                %                         plot(ROC(:,1),ROC_detrend,'b-')
                %                         axis([min_ROC_time max_ROC_time 0 1])
                %                         plot([0 0],[0 1],'k--')
                %                         plot([min_ROC_time max_ROC_time],[.5 .5],'k--')
                %                         plot([min_ROC_time max_ROC_time],[.6 .6],'k-.')
                %                         plot([min_ROC_time max_ROC_time],[.4 .4],'k-.')
                onset_latency_detrended = find_onset_latency_ROC(ROC_time,ROC_detrend,[30 100],[.6 8 10],0);
                if ~isnan(onset_latency)
                    xline(min(onset_latency, onset_latency_detrended), '-.k')
                    text(0,5,['ROC =',num2str(min(onset_latency, onset_latency_detrended))]);
                end
            elseif onset_latency ~= 0
                xline(onset_latency, '-.k')
                text(0,5,['ROC =',num2str(onset_latency)]);
            end
        end
        %         end
    end

end

end
