% p = r_diff_jo(R1,N1,R2,N2)
%   two-side significance level of difference between 
%   two correlation coeffients R1 and R2 based on 
%   N1 and N2 data points, respectively. (N1>=10, N2>=10)
%   
%   Jeroen Goossens

function p=r_diff_jo(R1,N1,R2,N2)

% significance level of difference between two 
% linear correlation coefficients R1 and R2.
% valid for moderate number of data points (N>=10)
% see: Press et. al.
%      Numerical Recipes 2ed edition p 636 

Z1 = 1/2.*log( (1+R1)./(1-R1) );  % Fisher's Z-transformation
Z2 = 1/2.*log( (1+R2)./(1-R2) );  

num = abs(Z1-Z2);
den = sqrt(2) .* sqrt( 1./(N1-3) + 1./(N2-3) );

p = erfc(num./den);



