function plot_SRT(SRT,plot_name,minRT, maxRT);
% simple plotting function for SRTs
% Also puts mean and stddev in hist plot
% minRT and maxRT are options. Imposes 80 and 600 by default
% Usage plot_SRT(SRT,plot_name,minRT,maxRT)

% Remove NaN and impose >80 < 600 limits

SRT = SRT(~isnan(SRT));
if nargin < 3; minRT = 80; maxRT = 600; end

SRT = SRT(SRT >= minRT & SRT <= maxRT);


EDGES = [0:10:maxRT];
N = histc(SRT,EDGES)/length(SRT);

hold on;
title(plot_name)
if size(N,1)>size(N,2); N = N'; end;
if ~isnan(SRT)
    bar(EDGES,N,'histc')
    axis([0 maxRT 0 0.3])
    text(25,0.25,strcat('Mean = ',num2str(mean(SRT)),' +/- ',num2str(std(SRT)), ', n = ',num2str(length(SRT))));
end