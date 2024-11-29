function [N,EDGES_out]=plot_histogram_raw(datain,EDGES,plot_name,stats_flag,color);
% simple plotting function for any data in a histogram program
% defined by EDGES input
% usage:
%   plot_histogram(datain,EDGES,plot_name,stats_flag)
% Stats_flag defaulted to 1 -- will show mean +/- std
% plot_name defaulted to '', unless provided
% EDGES typically [0:10:500]
% This plots raw numbers, not proportions

if nargin < 5
    color = 'w';
end
if nargin < 4
    stats_flag = 1;
end
if nargin < 3
    plot_name = '';
end

%datain = datain(~isnan(datain));
%datain = datain(outlier_jo(datain,3));

N = histc(datain,EDGES);
stepsize = EDGES(2)-EDGES(1);


hold on;
if nargin >= 3 & ~isempty(plot_name);
    title(plot_name)
end
bar(EDGES+stepsize/2,N,color)
axis([EDGES(1) EDGES(end) 0 max(N)*1.1])
start_text = EDGES(1) + 0.1 * (EDGES(end)-EDGES(1));
if stats_flag
    datain = datain(~isnan(datain));
    text(start_text,0.25,strcat('Mean = ',num2str(mean(datain)),' +/- ',num2str(std(datain)), ', n = ',num2str(length(datain))));
end

EDGES_out = EDGES+stepsize/2;
