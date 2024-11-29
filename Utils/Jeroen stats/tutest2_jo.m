% [p,t,df]=tutest2_jo(Mx,SDx,Nx,My,SDy,Ny)
%  Student's t-test on means Mx & My with standard deviations
%  SDx & SDy based on Nx and Ny observations.
%  The test indicates whether the data sets have significantly
%  different means. This procedure allows that the data sets X
%  and Y have unequal variances. Use the F-test to see whether 
%  this is the case. Small values of p indicate a significant 
%  difference. See also Press et al, pp 617
%
%   p   = signifficance level
%   t   = Student's t-value
%   df  = deegres of freedom
%
% Jeroen Goossens

function [p,t,df]=tutest2_jo(Mx,SDx,Nx,My,SDy,Ny)

VarX = SDx.^2;
VarY = SDy.^2;

% standard error
se = sqrt(VarX/Nx + VarY/Ny);

% t-value
t = (Mx-My)/se;

% degrees of freedom
df = (VarX/Nx + VarY/Ny)^2 / ( (VarX/Nx)^2/(Nx-1) + (VarY/Ny)^2/(Ny-1) );

% significance level
Xb = df/(df + t^2);  % see press et al, pp 228 eqn 6.4.9
Ab = df/2;
Bb = 1/2 ;
p  = betainc(Xb,Ab,Bb);

