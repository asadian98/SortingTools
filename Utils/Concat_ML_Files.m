function Z = Concat_ML_Files(ses_name, Physio_Info, ML_tail_list, RP_tail_list, animal_name, Folders, conditionCodeMap, strobecodesMap)

% Originally written by Dr. Aaron Cecala and extensively edited by Amirhossein Asadian
% To modify the target direction, adjust the corresponding function at the end of this script.

disp(' ------ Getting Z_ML .......')

Z_tr_offset = 0;
predur_min = Physio_Info.predur_min;
predur_max = Physio_Info.predur_max;
postdur = Physio_Info.postdur;
% Check how direction is calculated for Grover and Belle
% Check the codes for each condition

tic

% Butterworth Filter for Eye Traces (provided by Tyler Peel)
[A,B]=butter(3,0.0918); % Butterworth coefficients for G. 3rd order lowpass filter with fc = 0.0918 fs, with fs = 1000, fc = 45.9. Emulates the old Usui and Amidror 1982 values.

% What do you want to make movement direction based on coordinated
% eye-reach trials? Note: Eye-only and hand-only trials would be based
% on the movement direction relative to the primary effector in those
% trials.
% movedirection_based_on_IEP = 1; % based on initial eye position
% movedirection_based_on_IHP = 0; % based on initial hand position

for tail_idx = 1:length(ML_tail_list)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%              Check ML and Ripple recorded Statecodes          %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ------------------------------ Get Monkey Logic data ----------------
    cd(Folders.ML_Behav_Folder)
    ML_filename = [ses_name, ML_tail_list{tail_idx}, '.bhv2'];
    disp([ses_name, ML_tail_list{tail_idx}])
    [ML_behavioural_data,ML_Config,ML_TrialRecord] = mlread(ML_filename);
    st_ML = [];
    for i = 1:size(ML_behavioural_data, 2)
        st_ML = [st_ML; ML_behavioural_data(i).BehavioralCodes.CodeNumbers];
    end

    start_tr_num_list = 1;
    end_tr_num_list = size(ML_behavioural_data,2);

    % ------------------------------ Get Statecode ------------------------

    if(~isempty(RP_tail_list))
        statecodes_all = get_ML_statecodes(Folders, animal_name, ses_name, RP_tail_list, tail_idx);

        % Check if states codes are the same:
        if(statecodes_all(1) == statecodes_all(2))
            id = 1:2:length(statecodes_all);
            len = length(statecodes_all)/2;
            %         plot(st_ML); hold on; plot(statecodes_all(1:2:end))
        else
            id = 1:length(statecodes_all);
            len = length(statecodes_all);
            %         plot(st_ML); hold on; plot(statecodes_all)
        end

        if(length(st_ML) == len)
            disp('Same size')
            start_tr_num_list = 1;
            end_tr_num_list = size(ML_behavioural_data,2);
        else
            disp('Check ML size')

            [c, lags] = xcorr(st_ML, statecodes_all(id));
            [~, I] = max(abs(c));
            lag = lags(I);

            end_tr_num_list = size(ML_behavioural_data,2);

            start_tr_num_list = length(find(st_ML(1:lag) == 9)) + 1;
            Z_tr_offset = Z_tr_offset - length(find(st_ML(1:lag) == 9));

            disp(['The lag is: ', num2str(lag), ' means ', num2str(length(find(st_ML(1:lag) == 9))), ' trials']);

            cd(Folders.save_dir)
            % write a text file for errors
            fileID = fopen('Errors.txt', 'a');

            if(lag == 0)
                if length(st_ML) < len
                    fprintf(fileID, '%s\n', ['In ', ses_name, ML_tail_list{tail_idx}, ' strobe codes did not match, st_ML is shorter!']);
                else
                    st_ML_sh = st_ML(1:(len));
                    num9 = find(st_ML_sh == 9);
                    num18 = find(st_ML_sh == 18);
                    end_tr_num_list = end_tr_num_list - (length(num9) - length(num18) + length(find(st_ML(len + 1:end) == 9)));
                    fprintf(fileID, '%s\n', ['In ', ses_name, ML_tail_list{tail_idx}, ' strobe codes did not match, removed ', num2str(length(num9) - length(num18) + length(find(st_ML(len + 1:end) == 9))), ' trials from the end']);
                end
            elseif(lag < 0)
                fprintf(fileID, '%s\n', ['In ', ses_name, ML_tail_list{tail_idx}, ' strobe codes did not match --- Lag is nagative!!!!!, check Phy output with your data (ex. Analyze function)']);
            else
                fprintf(fileID, '%s\n', ['In ', ses_name, ML_tail_list{tail_idx}, ' strobe codes did not match, check Phy output with your data (ex. Analyze function)']);
            end

            fclose(fileID);

            cd(Folders.ML_Behav_Folder)

        end
    end
    Z.tr_sep(tail_idx) = end_tr_num_list - start_tr_num_list + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%               Create "Z" Structure              %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PARADIGM AND FILE SPECIFIC INFORMATION HERE

    parainfo.perfcodes{1,1} = 1;    parainfo.perfcodes{1,2} = 'Correct';
    parainfo.perfcodes{2,1} = 0;    parainfo.perfcodes{2,2} = 'Incorrect';
    parainfo.perfcodes{3,1} = -1;   parainfo.perfcodes{3,2} = 'Auto rejected due to RT outside min/max';
    parainfo.perfcodes{4,1} = -2;   parainfo.perfcodes{4,2} = 'Subject broke fixation';
    parainfo.perfcodes{5,1} = -99;  parainfo.perfcodes{5,2} = 'User rejected';

    keys = strobecodesMap.keys;
    for keys_idx = 1:length(keys)
        key = keys{keys_idx};
        value = strobecodesMap(key);  % Retrieve the struct
        parainfo.strobecodes{key,1} = value.code;  parainfo.strobecodes{value.code,2} = value.description;
    end

    start_gap_code = parainfo.strobecodes{12,1};
    end_gap_code = parainfo.strobecodes{13,1};

    ML_analoginfo.Eye{1,1} = 'Horizontal Eye'; ML_analoginfo.eye{1,2} = 'Vertical Eye';
    ML_analoginfo.Eye2{1,1} = 'EMPTY'; ML_analoginfo.eye2{1,2} = 'EMPTY';
    ML_analoginfo.EyeExtra{1,1} = 'EMPTY';  ML_analoginfo.EyeExtra{1,2} = 'EMPTY';
    ML_analoginfo.Joystick{1,1} = 'EMPTY'; ML_analoginfo.Joystick{1,2} = 'EMPTY';
    ML_analoginfo.Joystick2{1,1} = 'EMPTY'; ML_analoginfo.Joystick2{1,2} = 'EMPTY';
    ML_analoginfo.Touch{1,1} = 'Horizontal Hand'; ML_analoginfo.Touch{1,2} = 'Vertical Hand';
    ML_analoginfo.Mouse{1,1} = 'EMPTY'; ML_analoginfo.Mouse{1,2} = 'EMPTY';
    ML_analoginfo.KeyInput{1,1} = 'EMPTY'; ML_analoginfo.KeyInput{1,2} = 'EMPTY';
    ML_analoginfo.Photodiode{1,1} = 'EMPTY';
    ML_analoginfo.Photodiode{1,1} = 'EMPTY';

    ML_analoginfo.General.gen1{1,1}  = 'Photodiode associated with relevant target onset';
    ML_analoginfo.General.gen2{1,1}  = 'EMPTY';
    ML_analoginfo.General.gen3{1,1}  = 'EMPTY';
    ML_analoginfo.General.gen4{1,1}  = 'EMPTY';
    ML_analoginfo.General.gen5{1,1}  = 'EMPTY';
    ML_analoginfo.General.gen6{1,1}  = 'EMPTY';
    ML_analoginfo.General.gen7{1,1}  = 'EMPTY';
    ML_analoginfo.General.gen8{1,1}  = 'EMPTY';
    ML_analoginfo.General.gen9{1,1}  = 'EMPTY';
    ML_analoginfo.General.gen10{1,1} = 'EMPTY';

    ML_analoginfo.Button.btn1{1,1}  = 'EMPTY';
    ML_analoginfo.Button.btn2{1,1}  = 'EMPTY';
    ML_analoginfo.Button.btn3{1,1}  = 'EMPTY';
    ML_analoginfo.Button.btn4{1,1}  = 'EMPTY';
    ML_analoginfo.Button.btn5{1,1}  = 'EMPTY';
    ML_analoginfo.Button.btn6{1,1}  = 'EMPTY';
    ML_analoginfo.Button.btn7{1,1}  = 'EMPTY';
    ML_analoginfo.Button.btn8{1,1}  = 'EMPTY';
    ML_analoginfo.Button.btn9{1,1}  = 'EMPTY';
    ML_analoginfo.Button.btn10{1,1} = 'EMPTY';

    Z.parainfo = parainfo;
    Z.ML_analoginfo = ML_analoginfo;

    %% Specify and Extract Variables from MonkeyLogic File

    % Some files like 624C doesn't have the same number for Z.trials and
    % strobe code starts !!!
    %     if(tail_idx == 2)
    %         tr_num_list = 3:size(ML_behavioural_data,2);
    %         Z_tr_offset = Z_tr_offset - 2;
    %     else
    %         tr_num_list = 1:size(ML_behavioural_data,2);
    %     end

    tr_num_list = start_tr_num_list:end_tr_num_list;
    for tr_num = tr_num_list
        disp(['Processing MonkeyLogic Behavioural Data for trial number ',num2str(tr_num)])
        Z.TrialNumber(tr_num + Z_tr_offset) = ML_behavioural_data(tr_num).Trial;  % Actual Trial Number
        Z.Block(tr_num + Z_tr_offset) = ML_behavioural_data(tr_num).Block;        % Block Number
        Z.TrialWithinBlock(tr_num + Z_tr_offset) = ML_behavioural_data(tr_num).TrialWithinBlock; % Trial Within Block

        % Convert  TrialError in ML File to perfcodes defined above
        temp_TrialError = ML_behavioural_data(tr_num).TrialError;
        if temp_TrialError == 0 % Correct Trial
            Z.TrialError(tr_num + Z_tr_offset) = 1;
        elseif temp_TrialError == 1
            Z.TrialError(tr_num + Z_tr_offset) = 0;
        elseif temp_TrialError == 3
            Z.TrialError(tr_num + Z_tr_offset) = -2;
        elseif temp_TrialError == 4
            Z.TrialError(tr_num + Z_tr_offset) = -2;
        elseif temp_TrialError == 9
            Z.TrialError(tr_num + Z_tr_offset) = -2;
        else
            disp(['Trial error not classified for trial number: ',num2str(tr_num)])
        end
        Z.StrobeCodeNumbers{tr_num + Z_tr_offset} = ML_behavioural_data(tr_num).BehavioralCodes.CodeNumbers; % Strobe Code Times
        Z.StrobeCodeTimes{tr_num + Z_tr_offset} = ML_behavioural_data(tr_num).BehavioralCodes.CodeTimes; % Strobe Code Numbers

        % Get all Trial Variables intantiated within a given trial
        fields_UserVars = fieldnames(ML_behavioural_data(tr_num).UserVars);
        for fn_UV = 1:length(fields_UserVars)
            TrialVars{fn_UV,1} = fields_UserVars{fn_UV};
            TrialVars{fn_UV,2} = getfield(ML_behavioural_data(tr_num).UserVars,fields_UserVars{fn_UV});
        end
        Z.TrialVars{tr_num + Z_tr_offset} = TrialVars;
        % Get all object parameters instantiated in a give trial
        fields_SceneParameters = fieldnames(ML_behavioural_data(tr_num).ObjectStatusRecord.SceneParam);
        for fn_SP = 1:length(fields_SceneParameters)
            SceneParameters{fn_SP,1} = fields_SceneParameters{fn_SP};
            SceneParameters{fn_SP,2} = getfield(ML_behavioural_data(tr_num).ObjectStatusRecord.SceneParam,fields_SceneParameters{fn_SP});
        end
        Z.SceneParameters{tr_num + Z_tr_offset} = SceneParameters;
        % Catalog original ML Trial Type and convert ML Trial Type to
        % correct "parainfo.condition"
        % NOTE: This list of trial numbers may change with each
        % paradigm/experiment!!!
        if isfield(ML_behavioural_data(tr_num).UserVars,'Trial_Type')
            ML_TrialTypeNumber = getfield(ML_behavioural_data(tr_num).UserVars,'Trial_Type');
            Z.ML_TrialTypeNumber(tr_num + Z_tr_offset) = ML_TrialTypeNumber;

            % __________________________ Set Condition number -------------
            % You are giving the function a container map which includes
            % task code (the one you defined in ML), condition number (the
            % one you want to use in your analysis) and description
            key = Z.ML_TrialTypeNumber(tr_num + Z_tr_offset);
            if isKey(conditionCodeMap, key)
                value = conditionCodeMap(key);  % Retrieve the struct
                Z.condition(tr_num + Z_tr_offset) = value.condition;
            else
                disp('You have an undefined condition in your Z structure')
            end

            parainfo.paradigm = 'Plinko SC Recordings';
            parainfo.conditions{value.condition,1} = value.condition; parainfo.conditions{value.condition,2} =  value.description; % Saccade only task (trial type number: 502)

        end

        % Get all ML AnalogInputs within a given trial
        fields_AnalogData = fieldnames(ML_behavioural_data(tr_num).AnalogData);
        for fn_AD = 1:length(fields_AnalogData)
            AnalogData{fn_AD,1} = fields_AnalogData{fn_AD};
            AnalogData{fn_AD,2} = getfield(ML_behavioural_data(tr_num).AnalogData,fields_AnalogData{fn_AD});
        end
        Z.ML_AnalogData{tr_num + Z_tr_offset} = AnalogData;

        % Splice out the data you want
        % Timeframe to splice data
        Z.time_re_target1 = -predur_min:postdur;
        Z.time_re_target2 = -predur_max:postdur;

        % find first time that photodiode went high based on
        % thresholded value
        photothreshold = 0.5; %mV
        photodiode_data = Z.ML_AnalogData{1,tr_num + Z_tr_offset}{11,2}.Gen1;
        photodiode_data_greaterthanthreshold = find(photodiode_data(71:end) > photothreshold,1) + 70;
        if ~isempty(photodiode_data_greaterthanthreshold)
            Z.target_on_diode(tr_num + Z_tr_offset) = photodiode_data_greaterthanthreshold;
        else
            Z.target_on_diode(tr_num + Z_tr_offset) = NaN;
        end

        if(Z.condition(tr_num + Z_tr_offset) == 5 || Z.condition(tr_num + Z_tr_offset) == 6)
            predata_duration = predur_max;
        else
            predata_duration = predur_min;
        end

        % Get spliced Eye, Touch, and Photodiode

        if Z.TrialError(tr_num + Z_tr_offset) == 1
            Z.EyeX_Raw_Diode{tr_num + Z_tr_offset} = Z.ML_AnalogData{1,tr_num + Z_tr_offset}{2,2}(round(Z.target_on_diode(tr_num + Z_tr_offset))-predata_duration:end,1);
            Z.EyeY_Raw_Diode{tr_num + Z_tr_offset} = Z.ML_AnalogData{1,tr_num + Z_tr_offset}{2,2}(round(Z.target_on_diode(tr_num + Z_tr_offset))-predata_duration:end,2); % round(Z.target_on_diode(tr_num + Z_tr_offset))+postdata_duration

            % Filtered Eye
            Z.EyeX_Filt_Diode{tr_num + Z_tr_offset} = filtfilt(A,B,Z.EyeX_Raw_Diode{tr_num + Z_tr_offset});
            Z.EyeY_Filt_Diode{tr_num + Z_tr_offset} = filtfilt(A,B,Z.EyeY_Raw_Diode{tr_num + Z_tr_offset});

            % Photodiode
            Z.Photodiode_Diode{tr_num + Z_tr_offset} = photodiode_data(round(Z.target_on_diode(tr_num + Z_tr_offset))-predata_duration:end);

            % Touchscreen
            if (~isempty(Z.ML_AnalogData{1, 1}{7, 2}) && isempty(Z.ML_AnalogData{1, 1}{8, 2})) % touch data, no mouse data
                Z.TouchX_Diode{tr_num + Z_tr_offset} = Z.ML_AnalogData{1, tr_num + Z_tr_offset}{7, 2}(round(Z.target_on_diode(tr_num + Z_tr_offset))-predata_duration:end,1); % round(Z.target_on_diode(tr_num + Z_tr_offset))+postdata_duration
                Z.TouchY_Diode{tr_num + Z_tr_offset} = Z.ML_AnalogData{1, tr_num + Z_tr_offset}{7, 2}(round(Z.target_on_diode(tr_num + Z_tr_offset))-predata_duration:end,2); % round(Z.target_on_diode(tr_num + Z_tr_offset))+postdata_duration
            elseif (isempty(Z.ML_AnalogData{1, 1}{7, 2}) && ~isempty(Z.ML_AnalogData{1, 1}{8, 2})) % no touch data, mouse data
                Z.TouchX_Diode{tr_num + Z_tr_offset} = Z.ML_AnalogData{1, tr_num + Z_tr_offset}{8, 2}(round(Z.target_on_diode(tr_num + Z_tr_offset))-predata_duration:end,1); % round(Z.target_on_diode(tr_num + Z_tr_offset))+postdata_duration
                Z.TouchY_Diode{tr_num + Z_tr_offset} = Z.ML_AnalogData{1, tr_num + Z_tr_offset}{8, 2}(round(Z.target_on_diode(tr_num + Z_tr_offset))-predata_duration:end,2); % round(Z.target_on_diode(tr_num + Z_tr_offset))+postdata_duration
            end

            % Set GAP - TARGET LOCATION - TARGET_COLOUR
            Z = setTaskconstraint(Z, tr_num, Z_tr_offset, ML_behavioural_data, end_gap_code, start_gap_code, animal_name);

        else
            % Raw Eye
            Z.EyeX_Raw_Diode{tr_num + Z_tr_offset} = NaN(1,length(Z.time_re_target1))';
            Z.EyeY_Raw_Diode{tr_num + Z_tr_offset} = NaN(1,length(Z.time_re_target1))';

            % Filtered Eye
            Z.EyeX_Filt_Diode{tr_num + Z_tr_offset} = NaN(1,length(Z.time_re_target1))';
            Z.EyeY_Filt_Diode{tr_num + Z_tr_offset} = NaN(1,length(Z.time_re_target1))';

            % Photodiode
            Z.Photodiode_Diode{tr_num + Z_tr_offset} = NaN(1,length(Z.time_re_target1))';

            % Touchscreen
            Z.TouchX_Diode{tr_num + Z_tr_offset} = NaN(1,length(Z.time_re_target1))';
            Z.TouchY_Diode{tr_num + Z_tr_offset} = NaN(1,length(Z.time_re_target1))';

            % Reaction Time
            Z.Reach_RT_Diode(tr_num + Z_tr_offset) = NaN;

            % Target Characteristics
            Z.Target_Colour{tr_num + Z_tr_offset} = 'NONE';
            Z.Target_Location{1,tr_num + Z_tr_offset} = NaN;
            Z.Target_Location{2,tr_num + Z_tr_offset} = 'NONE';

            Z.GapDur(tr_num + Z_tr_offset) = nan;
        end

    end

    Z_tr_offset = Z_tr_offset + tr_num;

end

Z = detect_eye_arm_Ripple(Z, predur_min, predur_max);
Z.Physio_Info = Physio_Info;

disp(' ------ Saving Z_ML .......')
cd(Folders.save_dir)
Z_ML = Z;
save(['Z_', ses_name, '_ML'], 'Z_ML')

toc

    function Z = setTaskconstraint(Z, tr_num, Z_tr_offset, ML_behavioural_data, end_gap_code, start_gap_code, animal_name)

        if(isfield(ML_behavioural_data(tr_num).UserVars, 'Gap_Duration'))
            gap = ML_behavioural_data(tr_num).UserVars.Gap_Duration;

            if(gap < 50)
                gap = 0;
            elseif(75 < gap & gap < 125)
                gap = 100;
            elseif(125 < gap & gap < 175)
                gap = 150;
            elseif(175 < gap & gap < 225)
                gap = 200;
            elseif(225 < gap & gap < 275)
                gap = 250;
            elseif(275 < gap & gap < 325)
                gap = 300;
            elseif(gap > 300)
                gap = 300;
            end

            Z.GapDur(tr_num + Z_tr_offset) = gap;

        elseif(isfield(ML_behavioural_data(tr_num).UserVars, 'Gap_Dur'))

            gap = ML_behavioural_data(tr_num).UserVars.Gap_Dur;

            if(gap < 50)
                gap = 0;
            elseif(75 < gap & gap < 125)
                gap = 100;
            elseif(125 < gap & gap < 175)
                gap = 150;
            elseif(175 < gap & gap < 225)
                gap = 200;
            elseif(225 < gap & gap < 275)
                gap = 250;
            elseif(275 < gap & gap < 325)
                gap = 300;
            elseif(gap > 300)
                gap = 300;
            end

            Z.GapDur(tr_num + Z_tr_offset) = gap;
        else
            Z.GapDur(tr_num + Z_tr_offset) = nan;
        end

        % Get Target Location/Colour Information
        if Z.condition(tr_num + Z_tr_offset) == 1 || Z.condition(tr_num + Z_tr_offset) == 7 || Z.condition(tr_num + Z_tr_offset) == 8 %502

            if(isnan(Z.GapDur(tr_num + Z_tr_offset)))
                gap = ML_behavioural_data(tr_num).BehavioralCodes.CodeTimes(find(ML_behavioural_data(tr_num).BehavioralCodes.CodeNumbers == end_gap_code)) - ML_behavioural_data(tr_num).BehavioralCodes.CodeTimes(find(ML_behavioural_data(tr_num).BehavioralCodes.CodeNumbers == start_gap_code));

                if(gap < 50)
                    gap = 0;
                elseif(75 < gap & gap < 125)
                    gap = 100;
                elseif(125 < gap & gap < 175)
                    gap = 150;
                elseif(175 < gap & gap < 225)
                    gap = 200;
                elseif(225 < gap & gap < 275)
                    gap = 250;
                elseif(275 < gap & gap < 325)
                    gap = 300;
                elseif(gap > 300)
                    gap = 300;
                end

                Z.GapDur(tr_num + Z_tr_offset) = gap;
            end

            if(strcmp(animal_name, 'Belle'))
                T0_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T0_X')),2));
                T0_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T0_Y')),2));
                T1_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T1_X')),2));
                T1_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T1_Y')),2));
            else
                T0_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'T0_X')),2));
                T0_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'T0_Y')),2));
                T1_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'T1_X')),2));
                T1_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'T1_Y')),2));
            end

            if T1_X > T0_X
                Z.Target_Location{1,tr_num + Z_tr_offset} = 1;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'Right';
            elseif T1_X < T0_X
                Z.Target_Location{1,tr_num + Z_tr_offset} = -1;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'Left';
            else
                Z.Target_Location{1,tr_num + Z_tr_offset} = 0;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'VertOnly';
            end

            if T1_Y > T0_Y
                Z.Target_Location{3,tr_num + Z_tr_offset} = 1;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'Up';
            elseif T1_Y < T0_Y
                Z.Target_Location{3,tr_num + Z_tr_offset} = -1;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'Down';
            else
                Z.Target_Location{3,tr_num + Z_tr_offset} = 0;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'HoriOnly';
            end


            Z.Target_Colour{tr_num + Z_tr_offset} = 'GREEN';
        end

        if Z.condition(tr_num + Z_tr_offset) == 3 % -- just Grover

            Z.GapDur(tr_num + Z_tr_offset) = 0;

            T0_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T0_X')),2));
            T0_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T0_Y')),2));
            T1_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T1_X')),2));
            T1_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T1_Y')),2));

            if T1_X > T0_X
                Z.Target_Location{1,tr_num + Z_tr_offset} = 1;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'Right';
            elseif T1_X < T0_X
                Z.Target_Location{1,tr_num + Z_tr_offset} = -1;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'Left';
            else
                Z.Target_Location{1,tr_num + Z_tr_offset} = 0;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'VertOnly';
            end

            if T1_Y > T0_Y
                Z.Target_Location{3,tr_num + Z_tr_offset} = 1;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'Up';
            elseif T1_Y < T0_Y
                Z.Target_Location{3,tr_num + Z_tr_offset} = -1;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'Down';
            else
                Z.Target_Location{3,tr_num + Z_tr_offset} = 0;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'HoriOnly';
            end

            Z.Target_Colour{tr_num + Z_tr_offset} = 'GREEN';
        end

        if Z.condition(tr_num + Z_tr_offset) == 4 % Just Grover
            Z.GapDur(tr_num + Z_tr_offset) = 0;


            T0_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T0_X')),2));
            T0_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T0_Y')),2));
            T1_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYEHAND_T1_X')),2));
            T1_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYEHAND_T1_Y')),2));

            if T1_X > T0_X
                Z.Target_Location{1,tr_num + Z_tr_offset} = 1;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'Right';
            elseif T1_X < T0_X
                Z.Target_Location{1,tr_num + Z_tr_offset} = -1;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'Left';
            else
                Z.Target_Location{1,tr_num + Z_tr_offset} = 0;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'VertOnly';
            end

            if T1_Y > T0_Y
                Z.Target_Location{3,tr_num + Z_tr_offset} = 1;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'Up';
            elseif T1_Y < T0_Y
                Z.Target_Location{3,tr_num + Z_tr_offset} = -1;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'Down';
            else
                Z.Target_Location{3,tr_num + Z_tr_offset} = 0;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'HoriOnly';
            end


            Z.Target_Colour{tr_num + Z_tr_offset} = 'Blue';
        end

        if Z.condition(tr_num + Z_tr_offset) == 2 % Just Belle


            if(isnan(Z.GapDur(tr_num + Z_tr_offset)))
                gap = ML_behavioural_data(tr_num).BehavioralCodes.CodeTimes(find(ML_behavioural_data(tr_num).BehavioralCodes.CodeNumbers == end_gap_code)) - ML_behavioural_data(tr_num).BehavioralCodes.CodeTimes(find(ML_behavioural_data(tr_num).BehavioralCodes.CodeNumbers == start_gap_code));

                if(gap < 50)
                    gap = 0;
                elseif(75 < gap & gap < 125)
                    gap = 100;
                elseif(125 < gap & gap < 175)
                    gap = 150;
                elseif(175 < gap & gap < 225)
                    gap = 200;
                elseif(225 < gap & gap < 275)
                    gap = 250;
                elseif(275 < gap & gap < 325)
                    gap = 300;
                elseif(gap > 300)
                    gap = 300;
                end

                Z.GapDur(tr_num + Z_tr_offset) = gap;
            end


            %                         if movedirection_based_on_IEP == 1 && movedirection_based_on_IHP == 0
            %                             T0_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T0_X')),2));
            %                             T0_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYE_T0_Y')),2));
            %                             T1_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYEHAND_T1_X')),2));
            %                             T1_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYEHAND_T1_Y')),2));
            %                         elseif movedirection_based_on_IEP == 0 && movedirection_based_on_IHP == 1
            %                             T0_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'HAND_T0_X')),2));
            %                             T0_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'HAND_T0_Y')),2));
            %                             T1_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYEHAND_T1_X')),2));
            %                             T1_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'EYEHAND_T1_Y')),2));
            %                         else
            %                             disp('You forgot to choose whether the movement direction is based on IEP or IHP!')
            %                         end


            T0_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'HAND_T0_X')),2));
            T0_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'HAND_T0_Y')),2));
            T1_X = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'HAND_T1_X')),2));
            T1_Y = cell2mat(Z.TrialVars{1, tr_num + Z_tr_offset}(find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset},'HAND_T1_Y')),2));

            if T1_X > T0_X
                Z.Target_Location{1,tr_num + Z_tr_offset} = 1;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'Right';
            elseif T1_X < T0_X
                Z.Target_Location{1,tr_num + Z_tr_offset} = -1;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'Left';
            else
                Z.Target_Location{1,tr_num + Z_tr_offset} = 0;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'VertOnly';
            end

            if T1_Y > T0_Y
                Z.Target_Location{3,tr_num + Z_tr_offset} = 1;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'Up';
            elseif T1_Y < T0_Y
                Z.Target_Location{3,tr_num + Z_tr_offset} = -1;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'Down';
            else
                Z.Target_Location{3,tr_num + Z_tr_offset} = 0;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'HoriOnly';
            end

            Z.Target_Colour{tr_num + Z_tr_offset} = 'BLUE';
        end

        if Z.condition(tr_num + Z_tr_offset) == 5 || Z.condition(tr_num + Z_tr_offset) == 6 % ETP
            if Z.TrialVars{1,tr_num + Z_tr_offset}{find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset}, 'REACH_MT2_X')),2} == 25 || Z.TrialVars{1,tr_num + Z_tr_offset}{find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset}, 'REACH_MT2_X')),2} == 20
                Z.Target_Location{1,tr_num + Z_tr_offset} = 1;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'Right';
                Z.Target_Location{3,tr_num + Z_tr_offset} = 0;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'HoriOnly';

            end
            if Z.TrialVars{1,tr_num + Z_tr_offset}{find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset}, 'REACH_MT2_X')),2} == -5 || Z.TrialVars{1,tr_num + Z_tr_offset}{find(strcmp(Z.TrialVars{1, tr_num + Z_tr_offset}, 'REACH_MT2_X')),2} == -10
                Z.Target_Location{1,tr_num + Z_tr_offset} = -1;
                Z.Target_Location{2,tr_num + Z_tr_offset} = 'Left';
                Z.Target_Location{3,tr_num + Z_tr_offset} = 0;
                Z.Target_Location{4,tr_num + Z_tr_offset} = 'HoriOnly';
            end
        end

    end

    function statecodes_all = get_ML_statecodes(Folders, animal_name, ses_name, RP_tail_list, tail_idx)
        % _________________________________________________________________________
        %         For sessions that you splitted the nev file not the monkey logic
        %         file you need to change 1 to  1:length(ML_tail_list) on the next line, you should also comment some lines in the getEvents and
        %         Ripple Extract signals script!
        
        statecodes_all = [];

        for tail_idx_RP = 1
            disp('Getting Digital Statecodes for Behavioural Analysis')
            completeFilePath_Digital = [Folders.Ripple_Raw_Data_Folder, animal_name, '\', [ses_name, RP_tail_list{tail_idx}], '.nev']; % change tail_idx for multi segments recordings
            % Open the file and extract some basic information
            disp('opening Digital datafile to acquire statecodes')
            completeFilePath_Digital
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

            statecodes_all = [statecodes_all, statecodes];
        end

    end


end