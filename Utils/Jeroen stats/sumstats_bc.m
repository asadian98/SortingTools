function [data_mean, data_std, data_stderr, data_n] = sumstats_bc(datain, display_opt)
% Simple function created by bdc to summarize statistics of datain

if nargin == 1
    display_opt = 0;
end

datain = datain(~isnan(datain));

data_mean = mean(datain);
data_std = std(datain);
data_stderr = stderr_bc(datain);
data_n = length(datain);
data_min = min(datain);
data_max = max(datain);

if display_opt == 1;    % Echo text
    disp(strcat('Mean=',num2str(data_mean)));
    disp(strcat('Std=',num2str(data_std)));
    disp(strcat('Stderr=',num2str(data_stderr)));
    disp(strcat('n=',num2str(data_n)));
    disp(strcat('min =',num2str(data_min)));
    disp(strcat('max =',num2str(data_max)));
    
end