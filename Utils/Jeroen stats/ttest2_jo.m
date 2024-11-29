% [p,t,df]=ttest2_jo(Ho,M,SD,N)
%   student's t-test 
%
%   p   = signifficance level
%   t   = Student's t-value
%   df  = degrees of freedom
%
% Jeroen Goossens

function [p,t,df]=ttest2_jo(Ho,M,SD,N)

Var = SD.^2;

% standard error
se = sqrt(Var/N);

% t-value
t = (Ho-M)/se;

% degrees of freedom
df = N-1;

% significance level
Xb = df/(df + t^2);  % see press et al, pp 228 eqn 6.4.9
Ab = df/2;
Bb = 1/2 ;
p  = beta(Xb,Ab,Bb);
