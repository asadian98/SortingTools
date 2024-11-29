function plot_patchplot_noNANentries(xdata,ydata,plotcolor,plot_name);
% Produces patch plot, with mean +/- stderr
% USAGE: plot_patchplot(xdata,ydata,plotcolor,plot_name);
% color defaults to black if not provided
% Added feature is that it removes any NaN entries before calculating mean
% and STDERR

if nargin < 3;plotcolor = 'k';end
if nargin > 3
    if ~isempty(plot_name);
        title(plot_name)
    end
end

Ydata_noNAN_mean = [];
Ydata_noNAN_stderr = [];
for i = 1:size(ydata,2);  % For every column
    data = ydata(:,i);  % Grab column
    data = data(~isnan(data));  % Remove Nan Values
    Ydata_noNAN_mean(i) = mean(data);
    Ydata_noNAN_stderr(i) = stderr_bc(data);
end


h = patch([xdata fliparray(xdata)],[(Ydata_noNAN_mean + Ydata_noNAN_stderr ) fliparray(Ydata_noNAN_mean - Ydata_noNAN_stderr)],'k');
set(h,'EdgeColor',plotcolor)
set(h,'FaceColor',plotcolor)
return
