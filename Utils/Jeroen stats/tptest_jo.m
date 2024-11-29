%[p,t,df] = tptest_jo(X,Y)
%
% Student's t-test on two sets of paired data in vectors X
% and Y. The test indicates whether the difference between
% the means is significant. Small values of p indicate a
% significant difference.  See also Press et al, pp 618.
%
%   p   = signifficance level
%   t   = Student's t-value
%   df  = deegres of freedom
%
% Jeroen Goossens

function [p,t,df] = tptest_jo(X,Y)

% number of data points
Nx = length(X);
Ny = length(Y);
if Nx~=Ny,
  disp('Data in X and Y are not paired !');
  return;
end;
N = Nx;

% degrees of freedom
df = N-1;

% mean of distributions
Mx = mean(X);
My = mean(Y);

% variance and covariance of distributions
CovM  = cov(X,Y);
VarX  = CovM(1,1);
VarY  = CovM(2,2);
CovXY = CovM(1,2);

% standard error of the difference
se = sqrt( (VarX + VarY - 2*CovXY)/N );

% t-value
t = (Mx-My)/se;

% significance level
Xb = df/(df + t^2);  % see press et al, pp 228 eqn 6.4.9
Ab = df/2;
Bb = 1/2 ;
p  = betainc(Xb,Ab,Bb);

