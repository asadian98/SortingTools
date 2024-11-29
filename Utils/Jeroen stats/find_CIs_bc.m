function [low_CI high_CI] = find_CIs_bc(A,CI)

% Return the CI for the incoming dataset. Doing this across columns (each
% row is a different timepoint)
%
% Default values:
%   CI = 20%, ie, return 40 and 60
%   repeats = 500       % Number of times to repeat bootstrap

if nargin < 2
    CI = .20;
end

if CI == -1;
    % Indication that we're trying to do a non-parametric equivalent of
    % standard_error...
    CI=.68/sqrt(size(A,1));
    
end

low_CI_cut = ((1-CI)/2);
high_CI_cut = (1-(1-CI)/2);
    
numrows = size(A,1);
numcols = size(A,2);

low_CI = zeros(1,numcols);
high_CI = zeros(1,numcols);

low_CI_index = round(low_CI_cut*numrows);
high_CI_index = round(high_CI_cut * numrows);

for i = 1:numcols    
    A_sorted = sort(A(:,i));
    
    low_CI(i) = A_sorted(low_CI_index);
    high_CI(i) = A_sorted(high_CI_index);
end

return




