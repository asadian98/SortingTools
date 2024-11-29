function [low_CI,high_CI]=plot_patchplot_99(xdata,ydata,plotcolor,plot_name);
% Produces patch plot, covering 99% CI
% USAGE: plot_patchplot(xdata,ydata,plotcolor,plot_name);
% color defaults to black if not provided

if nargin < 3;plotcolor = 'k';end
if nargin > 3
    if ~isempty(plot_name);
        title(plot_name)
    end
end
high_CI = [];
low_CI = [];

for i = 1:size(ydata,2)
    low_CI(i) = mean(ydata(:,i)) + std(ydata(:,i)) * tinv(0.005,size(ydata,2)-1);
    high_CI(i) = mean(ydata(:,i)) + std(ydata(:,i)) * tinv(0.995,size(ydata,2)-1);    
end
patch([xdata fliparray(xdata)],[high_CI fliparray(low_CI)],plotcolor)

return
