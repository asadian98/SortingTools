% [r,p]=pearsn_jo(X,Y)
%   Pearson's (or linear) correlation coefficient
%
%   r = correlation coefficient
%   p = significance level
%
%   Jeroen Goossens

function [r,p]=pearsn_jo(X,Y)

r = corrcoef(X,Y);
r = r(1,2);
n = length(X);
p = r_signif_jo(r,n);
