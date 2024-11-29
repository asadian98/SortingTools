function [p] = chi2calc(x,df)
%CHI2CALC - Returns the probability of a chi-square value
%
%   function with V degrees of freedom at the values in X.
%   The chi-square density function with DF degrees of freedom,
%   is the same as a gamma density function with parameters V/2 and 2.
%
%   The size of P is the common size of X and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   This program calculates the chi-square significance values for given 
%   degrees of freedom and the tail probability (type I error rate) for 
%   given observed chi-square statistic and degree of freedom.

% Molecular Biology & Evolution Toolbox, (C) 2006
% Author: James J. Cai
% Email: jamescai@hku.hk
% Website: http://bioinformatics.org/mbetoolbox/
% Last revision: 5/28/2005

p=1-chi2cdf(x,df);

% find critical value
% critical = chi2inv(1-alpha, df);
% chi2inv(1-0.05,1)=3.8415