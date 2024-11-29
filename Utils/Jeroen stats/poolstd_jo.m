%
% POOLSTD_JO=(MEAN,STD,N)
%   POOLSTD calculates the pooled standard deviations of M
%   standard deviations STD(1:M) of means MEAN(1:M) which
%   are based on N(1:M) observations.
%
function ans=poolstd_jo(x,s,n);
dv  = mean(x)-x;
ans = sqrt( 1/(sum(n)-1) * sum( (n-1).*s.*s + n.*dv.*dv ) );

