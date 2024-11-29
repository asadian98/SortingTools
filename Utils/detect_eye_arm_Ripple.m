function Z = detect_eye_arm_Ripple(Z, predur_min, predur_max)

% Function created by bdc on Sep 30 2022 to run a rudimentary saccade
% detection algoritm on Data recorded by Aaron Cecala for NHPs working in
% Plinko condition. Adds Z.movearray and Z.eyeRT fields
%
% Usage detect_eye_Ripple('Z_G_220624_F')

% cd(['C:\Users\CorneilLab\Desktop\AA\Z_Files\', animal_name])

% disp('Loading. This may take a while')
% load(datain)

Z.EyeRT = [];   % NOTE THAT THIS OVERWRITES ALL PREVIOUS MARKS


for i = 1:size(Z.TrialNumber,2) % For every trial in the Z datastructure

    if Z.TrialError(i) > 0  % Only run saccade detection on "correct" trials
        Eh = Z.EyeX_Filt_Diode{i}; % Horizontal Eye Position
        Ev = Z.EyeY_Filt_Diode{i}; % Vertical Eye Position
        dEh = diff(Eh')*1000;
        dEv = diff(Ev')*1000;

        if(Z.condition(i) == 5 || Z.condition(i) == 6)
            DiodeTime = predur_max;
        else
            DiodeTime = predur_min;
        end
        output=sacdetector(dEh,dEv,Eh',Ev',[100 50 50],DiodeTime);
        if ~isempty(output)
            Z.movearray{i}=output;

            if 0     % For development, show saccade markings for each trial
                close all; figure;
                subplot(2,1,1); hold on; plot(Eh);plot(Ev,'r-'); ZZ = ylim;
                for x = 1:size(Z.movearray{i},1)
                    plot([Z.movearray{i}(x,1) Z.movearray{i}(x,1)],[ZZ(1) ZZ(2)],'g-')
                    plot([Z.movearray{i}(x,2) Z.movearray{i}(x,2)],[ZZ(1) ZZ(2)],'r-')
                end

                subplot(2,1,2); hold on; plot(dEh); plot(dEv,'r-'); ZZ = ylim;
                for x = 1:size(Z.movearray{i},1)
                    plot([Z.movearray{i}(x,1) Z.movearray{i}(x,1)],[ZZ(1) ZZ(2)],'g-')
                    plot([Z.movearray{i}(x,2) Z.movearray{i}(x,2)],[ZZ(1) ZZ(2)],'r-')
                end

            end

            AAA = output(:,3);   % RTs relative to target onset
            AAA = AAA(AAA>0); % RTs for eye movements occurring after target onset
            Z.EyeRT(i,1:length(AAA)) = AAA;
        end
    end


end
% Plot histograms of first saccades after T onset
A = Z.EyeRT(:,1); A= A(A>0); figure; plot_SRT(A,'',0,500);

% --------------------------------- Add arm-movement onset ----------------
ArmRT = zeros(1, length(Z.condition));
for tr = 1:length(Z.condition)
    if(Z.TrialError(tr) == 1)
        noTouch = find(isnan(Z.ML_AnalogData{tr}{7, 2}));
        if(isempty(noTouch))
            ArmRT(tr) = -999;
        else
            noTouchAfterTON_Idx = find(noTouch > Z.target_on_diode(tr));
            ArmRT(tr) = noTouch(noTouchAfterTON_Idx(1)) -  Z.target_on_diode(tr);
            if(ArmRT(tr) > 800)
                ArmRT(tr) = -999;
            end
        end
    else
        ArmRT(tr) = -999;
    end
end
Z.ArmRT = ArmRT;

% Previous code by Aaron for calc arm movement onset
% Strobe RT Calc
% Horizontal "velocity" values
%                     d_touch_h_GTThreshold_indicies_Strobe = find(d_touch_h_Strobe >= pos_threshold | d_touch_h_Strobe <= neg_threshold);
%                     d_touch_h_post_CueToGo_index_Strobe = find(d_touch_h_GTThreshold_indicies_Strobe>Target_Emergence_Time+min_RT_cutoff);
%
%                     % "NaN" values
%                     touch_h_NaN_Strobe = find(isnan(d_touch_h_Strobe));
%                     touch_h_NaN_post_CueToGo_indicies_Strobe = find(touch_h_NaN_Strobe>Target_Emergence_Time+min_RT_cutoff);
%
%                     % Determining Horizontal Movement Onset
%                     if ~isempty(touch_h_NaN_post_CueToGo_indicies_Strobe) & ~isempty(d_touch_h_post_CueToGo_index_Strobe)
%                         if touch_h_NaN_Strobe(touch_h_NaN_post_CueToGo_indicies_Strobe(1))< d_touch_h_GTThreshold_indicies_Strobe(d_touch_h_post_CueToGo_index_Strobe(1))
%                             h_movement_onset_Strobe = touch_h_NaN_Strobe(touch_h_NaN_post_CueToGo_indicies_Strobe(1));
%                         elseif touch_h_NaN_Strobe(touch_h_NaN_post_CueToGo_indicies_Strobe(1))> d_touch_h_GTThreshold_indicies_Strobe(d_touch_h_post_CueToGo_index_Strobe(1))
%                             h_movement_onset_Strobe = d_touch_h_GTThreshold_indicies_Strobe(d_touch_h_post_CueToGo_index_Strobe(1));
%                         end
%
%                     elseif isempty(touch_h_NaN_post_CueToGo_indicies_Strobe) & ~isempty(d_touch_h_post_CueToGo_index_Strobe)
%                             h_movement_onset_Strobe = d_touch_h_GTThreshold_indicies_Strobe(d_touch_h_post_CueToGo_index_Strobe(1));
%
%                     elseif isempty(d_touch_h_post_CueToGo_index_Strobe) & ~isempty(touch_h_NaN_post_CueToGo_indicies_Strobe)
%                             h_movement_onset_Strobe = touch_h_NaN_Strobe(touch_h_NaN_post_CueToGo_indicies_Strobe(1));
%
%                     elseif isempty(touch_h_NaN_post_CueToGo_indicies_Strobe) & isempty(d_touch_h_post_CueToGo_index_Strobe)
%                             h_movement_onset_Strobe = 5000;
%                     end
%
%                     % Vertical "velocity" values
%                     d_touch_v_GTThreshold_indicies_Strobe = find(d_touch_v_Strobe >= pos_threshold | d_touch_v_Strobe <= neg_threshold);
%                     d_touch_v_post_CueToGo_index_Strobe = find(d_touch_v_GTThreshold_indicies_Strobe>Target_Emergence_Time+min_RT_cutoff);
%                     % "NaN" values
%                     touch_v_NaN_Strobe = find(isnan(d_touch_v_Strobe));
%                     touch_v_NaN_post_CueToGo_indicies_Strobe = find(touch_v_NaN_Strobe>Target_Emergence_Time+min_RT_cutoff);
%
%                     % Determining Vertical Movement Onset
%                     if ~isempty(touch_v_NaN_post_CueToGo_indicies_Strobe) & ~isempty(d_touch_v_post_CueToGo_index_Strobe)
%                         if touch_v_NaN_Strobe(touch_v_NaN_post_CueToGo_indicies_Strobe(1))< d_touch_v_GTThreshold_indicies_Strobe(d_touch_v_post_CueToGo_index_Strobe(1))
%                             v_movement_onset_Strobe = touch_v_NaN_Strobe(touch_v_NaN_post_CueToGo_indicies_Strobe(1));
%                         elseif touch_h_NaN_Strobe(touch_v_NaN_post_CueToGo_indicies_Strobe(1))> d_touch_v_GTThreshold_indicies_Strobe(d_touch_v_post_CueToGo_index_Strobe(1))
%                             v_movement_onset_Strobe = d_touch_v_GTThreshold_indicies_Strobe(d_touch_v_post_CueToGo_index_Strobe(1));
%                         end
%                     elseif isempty(touch_v_NaN_post_CueToGo_indicies_Strobe) & ~isempty(d_touch_v_post_CueToGo_index_Strobe)
%                             v_movement_onset_Strobe = d_touch_v_GTThreshold_indicies_Strobe(d_touch_v_post_CueToGo_index_Strobe(1));
%                     elseif isempty(d_touch_v_post_CueToGo_index_Strobe) & ~isempty(touch_v_NaN_post_CueToGo_indicies_Strobe)
%                             v_movement_onset_Strobe = touch_v_NaN_Strobe(touch_v_NaN_post_CueToGo_indicies_Strobe(1));
%                     elseif isempty(touch_v_NaN_post_CueToGo_indicies_Strobe) & isempty(d_touch_v_post_CueToGo_index_Strobe)
%                             v_movement_onset_Strobe = 5000;
%                     end
%
%                     % Find Overall Reaction Time Value
%                     h_RT_Strobe = round(h_movement_onset_Strobe - Target_Emergence_Time);
%                     v_RT_Strobe = round(v_movement_onset_Strobe - Target_Emergence_Time);
%                         if (v_RT_Strobe <= h_RT_Strobe) & (v_RT_Strobe < max_RT_cutoff)
%                             O_RT_Strobe = v_RT_Strobe;
%                         elseif (v_movement_onset_Strobe >= h_RT_Strobe) & (h_RT_Strobe < max_RT_cutoff)
%                             O_RT_Strobe = h_RT_Strobe;
%                         else
%                             O_RT_Strobe = NaN;
%                         end
%
%                         Z.Reach_RT_Strobe(tr_num + Z_tr_offset) = O_RT_Strobe;

%             % Diode RT Calc
%             % Horizontal "velocity" values
%             d_touch_h_GTThreshold_indicies_Diode = find(d_touch_h_Diode >= pos_threshold | d_touch_h_Diode <= neg_threshold);
%             d_touch_h_post_CueToGo_index_Diode = find(d_touch_h_GTThreshold_indicies_Diode>Target_Emergence_Time+min_RT_cutoff);
%
%             % "NaN" values
%             touch_h_NaN_Diode = find(isnan(d_touch_h_Diode));
%             touch_h_NaN_post_CueToGo_indicies_Diode = find(touch_h_NaN_Diode>Target_Emergence_Time+min_RT_cutoff);
%
%             % Determining Horizontal Movement Onset
%             if ~isempty(touch_h_NaN_post_CueToGo_indicies_Diode) & ~isempty(d_touch_h_post_CueToGo_index_Diode)
%                 if touch_h_NaN_Diode(touch_h_NaN_post_CueToGo_indicies_Diode(1))< d_touch_h_GTThreshold_indicies_Diode(d_touch_h_post_CueToGo_index_Diode(1))
%                     h_movement_onset_Diode = touch_h_NaN_Diode(touch_h_NaN_post_CueToGo_indicies_Diode(1));
%                 elseif touch_h_NaN_Diode(touch_h_NaN_post_CueToGo_indicies_Diode(1))> d_touch_h_GTThreshold_indicies_Diode(d_touch_h_post_CueToGo_index_Diode(1))
%                     h_movement_onset_Diode = d_touch_h_GTThreshold_indicies_Diode(d_touch_h_post_CueToGo_index_Diode(1));
%                 end
%
%             elseif isempty(touch_h_NaN_post_CueToGo_indicies_Diode) & ~isempty(d_touch_h_post_CueToGo_index_Diode)
%                 h_movement_onset_Diode = d_touch_h_GTThreshold_indicies_Diode(d_touch_h_post_CueToGo_index_Diode(1));
%
%             elseif isempty(d_touch_h_post_CueToGo_index_Diode) & ~isempty(touch_h_NaN_post_CueToGo_indicies_Diode)
%                 h_movement_onset_Diode = touch_h_NaN_Diode(touch_h_NaN_post_CueToGo_indicies_Diode(1));
%
%             elseif isempty(touch_h_NaN_post_CueToGo_indicies_Diode) & isempty(d_touch_h_post_CueToGo_index_Diode)
%                 h_movement_onset_Diode = 5000;
%             end
%
%             % Vertical "velocity" values
%             d_touch_v_GTThreshold_indicies_Diode = find(d_touch_v_Diode >= pos_threshold | d_touch_v_Diode <= neg_threshold);
%             d_touch_v_post_CueToGo_index_Diode = find(d_touch_v_GTThreshold_indicies_Diode>Target_Emergence_Time+min_RT_cutoff);
%             % "NaN" values
%             touch_v_NaN_Diode = find(isnan(d_touch_v_Diode));
%             touch_v_NaN_post_CueToGo_indicies_Diode = find(touch_v_NaN_Diode>Target_Emergence_Time+min_RT_cutoff);
%
%             % Determining Vertical Movement Onset
%             if ~isempty(touch_v_NaN_post_CueToGo_indicies_Diode) & ~isempty(d_touch_v_post_CueToGo_index_Diode)
%                 if touch_v_NaN_Diode(touch_v_NaN_post_CueToGo_indicies_Diode(1))< d_touch_v_GTThreshold_indicies_Diode(d_touch_v_post_CueToGo_index_Diode(1))
%                     v_movement_onset_Diode = touch_v_NaN_Diode(touch_v_NaN_post_CueToGo_indicies_Diode(1));
%                 elseif touch_h_NaN_Diode(touch_v_NaN_post_CueToGo_indicies_Diode(1))> d_touch_v_GTThreshold_indicies_Diode(d_touch_v_post_CueToGo_index_Diode(1))
%                     v_movement_onset_Diode = d_touch_v_GTThreshold_indicies_Diode(d_touch_v_post_CueToGo_index_Diode(1));
%                 end
%             elseif isempty(touch_v_NaN_post_CueToGo_indicies_Diode) & ~isempty(d_touch_v_post_CueToGo_index_Diode)
%                 v_movement_onset_Diode = d_touch_v_GTThreshold_indicies_Diode(d_touch_v_post_CueToGo_index_Diode(1));
%             elseif isempty(d_touch_v_post_CueToGo_index_Diode) & ~isempty(touch_v_NaN_post_CueToGo_indicies_Diode)
%                 v_movement_onset_Diode = touch_v_NaN_Diode(touch_v_NaN_post_CueToGo_indicies_Diode(1));
%             elseif isempty(touch_v_NaN_post_CueToGo_indicies_Diode) & isempty(d_touch_v_post_CueToGo_index_Diode)
%                 v_movement_onset_Diode = 5000;
%             end
%
%             % Find Overall Reaction Time Value
%             h_RT_Diode = round(h_movement_onset_Diode - Target_Emergence_Time);
%             v_RT_Diode = round(v_movement_onset_Diode - Target_Emergence_Time);
%             if (v_RT_Diode <= h_RT_Diode) & (v_RT_Diode < max_RT_cutoff)
%                 O_RT_Diode = v_RT_Diode;
%             elseif (v_movement_onset_Diode >= h_RT_Diode) & (h_RT_Diode < max_RT_cutoff)
%                 O_RT_Diode = h_RT_Diode;
%             else
%                 O_RT_Diode = NaN;
%             end
%
%             Z.Reach_RT_Diode(tr_num + Z_tr_offset) = O_RT_Diode;


% UpFront Stuff
% Data before and after alignPt
%     predata_duration = 500;
%     postdata_duration = 500;
% Threshold for marking Touchscreen Data values
% pos_threshold = .4; % ok for Grover, but may need to modify this for other NHPs
% neg_threshold = -1*(pos_threshold);
% % Touchscreen Reaction Time min/max (ms)
% min_RT_cutoff = 50;
% max_RT_cutoff = 700;


end
% return


