function plot_patchplot_median(xdata,ydata,plotcolor,plot_name);
% Produces patch plot, with median +/- stderr
% USAGE: plot_patchplot(xdata,ydata,plotcolor,plot_name);
% color defaults to black if not provided

if nargin < 3;plotcolor = 'k';end
if nargin > 3
    if ~isempty(plot_name);
        title(plot_name)
    end
end
h = patch([xdata fliparray(xdata)],[(median(ydata) + stderr_bc(ydata) ) fliparray(median(ydata) - stderr_bc(ydata) )],'k');
set(h,'EdgeColor',plotcolor)
set(h,'FaceColor',plotcolor)
return
