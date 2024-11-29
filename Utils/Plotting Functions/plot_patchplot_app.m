function plot_patchplot_app(fig_handle, xdata,ydata,plotcolor,plot_name);
% Produces patch plot, with mean +/- stderr
% USAGE: plot_patchplot(xdata,ydata,plotcolor,plot_name);
% color defaults to black if not provided

if nargin < 3;plotcolor = 'k';end
if nargin > 4
    if ~isempty(plot_name)
        title(fig_handle, plot_name)
    end
end
h = patch(fig_handle, [xdata fliparray(xdata)],[(mean(ydata) + stderr_bc(ydata) ) fliparray(mean(ydata) - stderr_bc(ydata) )],'k');
set(h,'EdgeColor',plotcolor)
set(h,'FaceColor',plotcolor)
return
