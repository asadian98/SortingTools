function stderr = nanstderr_bd(data)
% Simple function to return standard error of 1D array

data = data(~isnan(data));

stderr = std(data)/sqrt(length(data));