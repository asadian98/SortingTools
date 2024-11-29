function R = plot_linreg(xdata,ydata,symbol,labels,axis_vals)
% Plotting function designed by BDC on Jan 17 2019 to do a simple linear
% regression and return the results 
% symbol, labels and axes are optional

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
    axis_vals(1) = min(xdata) - .1*min(xdata);
    axis_vals(2) = max(xdata) + .1*max(xdata);
    axis_vals(3) = min(ydata) - .1*min(ydata);
    axis_vals(4) = max(ydata) + .1*max(ydata);
end

hold on; 
plot(xdata,ydata,symbol);
R=linfit_jo(xdata,ydata);
pval_round = round(R(6)*1000)/1000;
rval_round = round(R(5)*100)/100;
if R(6) < 0.05
    lsline;
end
axis([axis_vals(1) axis_vals(2) axis_vals(3) axis_vals(4)])
axis square
T = text(axis_vals(1)*1.1,axis_vals(4)*.9,strcat('r=',num2str(rval_round),'p=',num2str(pval_round)));
color = symbol(1);
set(T,'color',color);
xlabel(xaxislabel)
ylabel(yaxislabel)
