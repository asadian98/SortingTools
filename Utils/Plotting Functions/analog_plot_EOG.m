function analog_plot_EOG(Ttime,datain,plottitle,scale,marks1,marks2);
% Generic plotting function for any analog data.
%   Ttime: time scale. Must match the columns of Data
%   Data: each row is a different trial
%   title: title of subplot
%   scale: min and max range of dataset
%   marks1: optional set of marks, relative to 0 of Ttime
%   marks2: optional set of marks, relative to 0 of Ttime

min_data = scale(1);
max_data = scale(2);
hold on;
title(plottitle)
for i = 1:size(datain,1)
    plot(Ttime,datain(i,:),'k-')
    if nargin == 6
        % PLOT MARKS 1 and MARKS 2
        % PLOT MARKS 1
        if ~isnan(marks1(i))
            plot(marks1(i),0,'r.');
        end
        if ~isnan(marks2(i))
            plot(marks2(i),0,'g.');
        end
    elseif nargin == 5
        % PLOT MARKS 1
        if ~isnan(marks1(i))
            plot(marks1(i),0,'r.');
        end
    end
end
line([0 0],[min_data max_data]);
axis([min(Ttime) max(Ttime) min_data max_data])
