function plot_CDF(datain,plotcolor,plottitle)
% Simple function designed by BDC on Oct 27 2022 to plot cumulative
% distribution function for datain

if nargin < 2
    plotcolor = 'b'
end
if nargin < 3
    plottitle = ''
end

hold on; 
title(plottitle)

datain_sorted = sort(datain);

for i = 1:length(datain_sorted)
    datain_sorted(i,2) = i/size(datain_sorted,1);
end

plot(datain_sorted(:,1),datain_sorted(:,2),'Color',plotcolor)


return

