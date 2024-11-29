function plot_patchplot_inline(axis, xdata,ydata,plotcolor,plot_name, facealpha)
% Produces patch plot, with mean +/- stderr
% USAGE: plot_patchplot(xdata,ydata,plotcolor,plot_name);
% color defaults to black if not provided

if nargin < 3;plotcolor = 'k';end
if nargin > 3
    if ~isempty(plot_name);
        title(axis, plot_name)
    end
end

if nargin < 5
    facealpha = 1;
end

h = patch(axis, [xdata fliparray(xdata)],[(nanmean(ydata) + stderr_bc(ydata)) fliparray(nanmean(ydata) - stderr_bc(ydata))],'k', 'FaceAlpha', facealpha, 'LineStyle',"none");
set(h,'EdgeColor',plotcolor)
set(h,'FaceColor',plotcolor)
return
