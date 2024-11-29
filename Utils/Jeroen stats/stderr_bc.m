function stderr = stderr_bd(data)
% Simple function to return standard error of 1D array
stderr = nanstd(data)/sqrt(length(data));