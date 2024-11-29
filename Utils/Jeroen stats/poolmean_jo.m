%
% POOLMEAN_JO=(MEAN,N)
%    POOLMEAN calculates the grand mean of M means
%    MEAN(1:M) which are based on N(1:M) observations.
%

function ans=poolmean_jo(x,n);
ans = 1/sum(n) * sum(n.*x);