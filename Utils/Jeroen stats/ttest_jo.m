%[p,t,df] = ttest_jo(X,Y)
%
% Student's t-test on two sets of data in vectors X and Y.
% The test indicates whether the data sets have significantly
% different means. It is assumed that the variances of data 
% sets X and Y are the same. Use the F-test to check if this 
% is true. Small values of p indicate a significant difference.
% See also Press et al, pp 616
%
%   p   = signifficance level
%   t   = Student's t-value
%   df  = deegres of freedom
%
% Jeroen Goossens

function [p,t,df] = ttest_jo(X,Y)

if ftest(X,Y)<0.05,
  disp('  Warning: Variances of X and Y may be different. ');
  disp('           You may use tutest as an alternative.  ');
end;

% number of data points
Nx = length(X);
Ny = length(Y);

% degrees of freedom
df = Nx + Ny -2 ;

% mean of distributions
Mx = mean(X);
My = mean(Y);

% variance of distributions
VarX = std(X)^2;
VarY = std(Y)^2;

% pooled variance
svar = ( (Nx-1)*VarX + (Ny-1)*VarY )/df;

% standard error of the difference
se = sqrt(svar*(1/Nx + 1/Ny));

% t-value
t = (Mx - My) / se;

% significance level
Xb = df/(df + t^2);  % see press et al, pp 228 eqn 6.4.9
Ab = df/2;
Bb = 1/2 ;
p  = betainc(Xb,Ab,Bb);

