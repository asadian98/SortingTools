function plot_double_histogram(datain,datain2,EDGES,plot_name,stats_flag,color);
% simple plotting function for any data in a histogram program
% defined by EDGES input
% datain2 plotted as inverse
% data plotted as frequency, not as observations
% stats_flag defaulted to turn on. Set to zero to suppress
% Color defaulted to white
% usage:
%   plot_double_histogram(datain,datain2,EDGES,plot_name,stats_flag,'w')


if nargin <= 5
    stats_flag = 1;
    color = 'w';
end
if nargin < 4
    plot_name = '';
end

stepsize = EDGES(2)-EDGES(1);

datain = datain(~isnan(datain));
datain = datain(outlier_jo(datain,3));
datain2 = datain2(~isnan(datain2));
datain2 = datain2(outlier_jo(datain2,3));

N = histc(datain,EDGES)/length(datain);
N2 = histc(datain2,EDGES)/length(datain2)*-1;

hold on;
if nargin >= 4 & ~isempty(plot_name);
    title(plot_name)
end
if isempty(N);
    ymax=.1;
else
    bar(EDGES+stepsize/2,N,color)
    ymax = max(N) * 1.1;
end
if isempty(N2);
    ymin = -.1;
else
    bar(EDGES+stepsize/2,N2,color)
    ymin = min(N2) * 1.1;
end
axis([EDGES(1) EDGES(end) ymin ymax])
start_text = EDGES(1) + 0.1 * (EDGES(end)-EDGES(1));
if stats_flag
    if ~isempty(N); text(start_text,ymax*.75,strcat('Mean = ',num2str(mean(datain)),' +/- ',num2str(std(datain)), ', n = ',num2str(length(datain)))); end
    if ~isempty(N2);text(start_text,ymin*.75,strcat('Mean = ',num2str(mean(datain2)),' +/- ',num2str(std(datain2)), ', n = ',num2str(length(datain2))));end
end
