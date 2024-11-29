function plot_scatterplot(xdata,ydata,symbol,labels,axis_vals)
% Plotting function designed by BDC on Nov 5 2019 to do a simple scatter
% plot and return results of paired t-test.
% The assumption is that the rows for xdata and ydata come for the same
% list of individuals, hence the ttest 
% symbol,axes are optional
% Usage plot_scatterplot(xdata,ydata,'ks',{'xlabel';'ylabel'},[0 100 0 100]);

if nargin < 3
    symbol = 'ks'
end
if nargin < 4
    xaxislabel = '';
    yaxislabel = '';
else
    xaxislabel = labels{1};
    yaxislabel = labels{2};
end
if nargin < 5
    axis_vals(1) = min([min(xdata) min(ydata)])-abs(.1*min([min(xdata) min(ydata)]));
    axis_vals(2) = max([max(xdata) max(ydata)])+.1*max([max(xdata) max(ydata)]);
    axis_vals(3) = min([min(xdata) min(ydata)])-abs(.1*min([min(xdata) min(ydata)]));
    axis_vals(4) = max([max(xdata) max(ydata)])+.1*max([max(xdata) max(ydata)]);
end

% Remove any rows with NAN values
X = [xdata ydata];
X = X(find(~isnan(X(:,1))&(~isnan(X(:,2)))),:);
xdata = X(:,1);ydata = X(:,2);

hold on;
plot(xdata,ydata,symbol)
axis square
axis([axis_vals(1) axis_vals(2) axis_vals(3) axis_vals(4)])
line([axis_vals(1) axis_vals(2)],[axis_vals(3) axis_vals(4)])
text(axis_vals(1)*1.1,axis_vals(2)*.9,num2str(tptest_jo(xdata,ydata)))
xlabel(xaxislabel)
ylabel(yaxislabel)

