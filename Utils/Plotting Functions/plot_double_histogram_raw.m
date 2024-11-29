function plot_double_histogram_raw(datain,datain2,EDGES,plot_name,stats_flag,color);
% simple plotting function for any data in a histogram program
% defined by EDGES input
% datain2 plotted as inverse
% data plotted as observations, not frequency
% stats_flag defaulted to turn on. Set to zero to suppress
% Color defaulted to white
% usage:
%   plot_double_histogram(datain,datain2,EDGES,plot_name,stats_flag,'w')


if nargin < 5
    stats_flag = 1;
    color = 'w';
end
if nargin < 4
    plot_name = '';
end

if nargin >= 4 & ~isempty(plot_name);
    title(plot_name)
end
stepsize = EDGES(2)-EDGES(1);

% Remove NANs and outliers
hold on;
if length(datain)>0
    datain = datain(~isnan(datain));
    datain = datain(outlier_jo(datain,3));
    N = histc(datain,EDGES);
    bar(EDGES+stepsize/2,N,color)
    start_text = EDGES(1) + 0.1 * (EDGES(end)-EDGES(1));
    if stats_flag
        text(start_text,max(N)*.75,strcat('Mean = ',num2str(mean(datain)),' +/- ',num2str(std(datain)), ', n = ',num2str(length(datain))));
    end
    ymax = max(N)*1.1;
else
    ymax = 1;
end

if length(datain2)>0    
    datain2 = datain2(~isnan(datain2));
    datain2 = datain2(outlier_jo(datain2,3));
    N2 = histc(datain2,EDGES)*-1;
    bar(EDGES+stepsize/2,N2,color)
    start_text = EDGES(1) + 0.1 * (EDGES(end)-EDGES(1));
    if stats_flag
        text(start_text,min(N2)*.75,strcat('Mean = ',num2str(mean(datain2)),' +/- ',num2str(std(datain2)), ', n = ',num2str(length(datain2))));
    end
    ymin = min(N2) * 1.1;
else
    ymin = -1;
end
axis([EDGES(1) EDGES(end) ymin ymax])
