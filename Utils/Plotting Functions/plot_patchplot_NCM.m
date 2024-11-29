function plot_patchplot_NCM(xdata,ydata,plotcolor,plot_name, facealpha)
% Produces patch plot, with mean +/- stderr
% USAGE: plot_patchplot(xdata,ydata,plotcolor,plot_name);
% color defaults to black if not provided

if nargin < 3;plotcolor = 'k';end
if nargin > 3
    if ~isempty(plot_name);
        title(plot_name)
    end
end
h = patch([xdata fliparray(xdata)],[(nanmean(ydata) + 2*stderr_bc(ydata)) fliparray(nanmean(ydata) - 2*stderr_bc(ydata))],'k', 'FaceAlpha', facealpha, 'LineStyle',"none");
set(h,'EdgeColor',plotcolor)
set(h,'FaceColor',plotcolor)
hold on

% plot(xdata, mean(ydata), 'color', 'k', 'linewidth', 1);
% set(h2,'EdgeColor',plotcolor)
% set(h2,'color',plotcolor)

return
