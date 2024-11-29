function plot_patchplot_noNAN(xdata,ydata,plotcolor,plot_name);
% Produces patch plot, with mean +/- stderr
% USAGE: plot_patchplot(xdata,ydata,plotcolor,plot_name);
% color defaults to black if not provided
% Added feature is that it removes any NaN rows first

if nargin < 3;plotcolor = 'k';end
if nargin > 3
    if ~isempty(plot_name);
        title(plot_name)
    end
end

ydata = ydata(~isnan(ydata(:,1)),:);

h = patch([xdata fliparray(xdata)],[(mean(ydata) + stderr_bc(ydata) ) fliparray(mean(ydata) - stderr_bc(ydata) )],'k');
set(h,'EdgeColor',plotcolor)
set(h,'FaceColor',plotcolor)
return
