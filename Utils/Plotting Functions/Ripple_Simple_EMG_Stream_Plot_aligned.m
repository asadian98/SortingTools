close all

%MODIFIED FROM AJC FUNCTION TO ALSO PLOT ONSET-ALIGNED AVERAGED

% Which Electrodes do you want?
Electrode_1 = 17; % should be an odd number and less than or equal to 19
Electrode_2 = Electrode_1+1;

% Scaling Factors for y-axis
Plot_y_upper = 1000;
Plot_y_lower = 1000;
monopolar_scale = 4;

% Total Time 
Total_Time_Minutes = round((length(ArrayData_Stream{Electrode_1,1})/2000)/60);
Total_Time_Seconds = round((length(ArrayData_Stream{Electrode_1,1})/2000));

% % Plot Time
% Number_of_Seconds_to_Plot = Total_Time_Seconds;
% Plotting_Convert_Seconds_to_Data_Points = Number_of_Seconds_to_Plot*2000;
% Starting_Point_of_Plot_in_Seconds = 1;
% Starting_Convert_Seconds_to_Data_Points = Starting_Point_of_Plot_in_Seconds*2000;

% Muscle Labels
if Electrode_1 == 1
    muscle_name = 'Pec Clav Head';
elseif Electrode_1== 3
    muscle_name = 'Deltoid Anterior';
elseif Electrode_1 == 5
    muscle_name = 'Deltoid Middle';
elseif Electrode_1 == 7
    muscle_name = 'Deltoid Posterior';
elseif Electrode_1 == 9
    muscle_name = 'Left Rectus';
elseif Electrode_1 == 11
    muscle_name = 'Right Rectus';
elseif Electrode_1 == 13
    muscle_name = 'Right OCi';
elseif Electrode_1 == 15
    muscle_name = 'Left OCi';
elseif Electrode_1 == 17
    muscle_name = 'Left Splenius';
elseif Electrode_1 == 19
    muscle_name = 'Right Splenius';
end
    
% Plot Electrode #1
subplot(3,1,1)
plot(ArrayData_Stream{Electrode_1,1});
ylim([-Plot_y_lower*monopolar_scale Plot_y_upper*monopolar_scale])
title(['Electrode',' ',num2str(Electrode_1)])
%xticks([0 (length(ArrayData_Stream{Electrode_1,1})/2) length(ArrayData_Stream{Electrode_1,1})])
%xticklabels({'0',num2str(Total_Time_Minutes/2) num2str(Total_Time_Minutes)})

% Plot Electrode #2
subplot(3,1,2)
plot(ArrayData_Stream{Electrode_2,1});
ylim([-Plot_y_lower*monopolar_scale Plot_y_upper*monopolar_scale])
ylabel('EMG (MicroV)')
title(['Electrode',' ',num2str(Electrode_2)])
%xticks([0 (length(ArrayData_Stream{Electrode_1,1})/2) length(ArrayData_Stream{Electrode_1,1})])
%xticklabels({'0',num2str(Total_Time_Minutes/2) num2str(Total_Time_Minutes)})

% Plot Differential Electrode
subplot(3,1,3)
plot(ArrayData_Stream{Electrode_1,1}-ArrayData_Stream{Electrode_2,1});
ylim([-Plot_y_lower Plot_y_upper])
title(['Differential Electrode: ', muscle_name])
xlabel('Time (minutes)')
%xticks([0 (length(ArrayData_Stream{Electrode_1,1})/2) length(ArrayData_Stream{Electrode_1,1})])
%xticklabels({'0',num2str(Total_Time_Minutes/2) num2str(Total_Time_Minutes)})


EventStream = ArrayData_Stream{end};

% CODE FOR DETECTING THE ONSETS
threshold = 1000;
A = find(EventStream > threshold);
B = A(diff(A) > 5);%Yields indices above threshold where the differences between adjacent points is greater than five. One such point per transition
crossings = zeros(1,length(B));
tic
crossings(1) = max(find(EventStream(1:B(1))<threshold));
for i = 2:length(B)
    i;
    crossings(i) = max(find(EventStream(crossings(i-1):B(i))<threshold))+crossings(i-1);
end
toc

% Use these crossings to calculate the event-aligned averages
pre = -100; post = 200;
time_re_event = pre:0.5:post;
E1 = [];
E2 = [];
DeltaE = [];

tic
for i = 1:length(crossings)
    time_pre = 2*crossings(i)+(pre*2);
    time_post = 2*crossings(i)+(post*2);
    E1(end+1,:) = abs(ArrayData_Stream{Electrode_1,1}(time_pre:time_post));
    E2(end+1,:) = abs(ArrayData_Stream{Electrode_1,1}(time_pre:time_post));
    DeltaE(end+1,:) = abs(ArrayData_Stream{Electrode_1,1}(time_pre:time_post)-ArrayData_Stream{Electrode_2,1}(time_pre:time_post));
end
toc

figure;
subplot(3,1,1); hold on;
plot_patchplot(time_re_event,E1)
subplot(3,1,2); hold on;
plot_patchplot(time_re_event,E2)
subplot(3,1,3); hold on;
plot_patchplot(time_re_event,DeltaE)
