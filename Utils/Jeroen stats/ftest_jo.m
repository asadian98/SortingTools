%[p,f,df1,df2]=ftest_jo(X,Y)
%  F-test for significanly different variances.
%  Small values of p indicate a significant difference.
%  See Press et. al. pp 619
%
%  p   = signifficance level
%  f   = f-value (ratio of variances)
%  df1 = degrees of freedom
%  df2 = degrees of freedom
%
%  Jeroen Goossens

function [p,f,df1,df2]=ftest_jo(X,Y)

% number of data points
Nx = length(X);
Ny = length(Y);

% mean of distributions
Mx = mean(X);
My = mean(Y);

% variance of distributions
VarX = std(X)^2;
VarY = std(Y)^2;

% Make F the ratio of the larger to the smaller one
if VarX>VarY,     
  f   = VarX/VarY;
  df1 = Nx-1;
  df2 = Ny-1;
else
  f   = VarY/VarX;
  df1 = Ny-1;
  df2 = Nx-1;
end;

% significance level
Xb = df2/(df2 + df1*f);  % see press et al, pp 228 eqn 6.4.9
Ab = df2/2; 
Bb = df1/2 ;
%p  = 2*beta(Xb,Ab,Bb);
p  = 2*betainc(Xb,Ab,Bb);
if (p>1), p = 2-p; end;