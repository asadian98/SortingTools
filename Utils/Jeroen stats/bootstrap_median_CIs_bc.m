function [low_CI high_CI] = bootstrap_median_CIs_bc(A,CI,repeats)

% Simple Monte-Carlo based bootstrapping procedure to estimate the
% confidence intervals of the median of the incoming distribution. 
%
% Default values:
%   CI = 50%, ie, return 25 and 75
%   repeats = 500       % Number of times to repeat bootstrap

if nargin < 2
    CI = .50;
    repeats = 500;
elseif nargin < 3
    repeats = 500;
end

low_CI = zeros(1,size(A,2));
high_CI = zeros(1,size(A,2));

for x = 1:size(A,2);
    
    bootstrap_median = zeros(1,repeats);
    
    distribution = A(:,x);
    
%    tic
    for i = 1:repeats
        
        bootstrap_median(i) = median(distribution(ceil(length(distribution) * rand(1,length(distribution)))));
        
    end
%    toc
    
    % Figure out CIs here
    high_CI_index = round((1-(1-CI)/2) * repeats);
    low_CI_index = round(((1-CI)/2) * repeats);
    bootstrap_median_sort = sort(bootstrap_median);
    
    high_CI(x) = bootstrap_median_sort(high_CI_index);
    low_CI(x) = bootstrap_median_sort(low_CI_index);
    
end

return




