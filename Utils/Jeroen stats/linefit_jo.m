% [Par,Std]=linefit_jo(X,Y)
%   fit a straight line through data set (X,Y)
%   and return the slope and intercept with their
%   standard deviations (using bootstrap method).
%
%   Par=[gain,offset],  Std=[std gain, std offset]
%
%   Jeroen Goossens

function [ParMean, ParStd]=linefit_jo(X,Y)

Par = bootstrp_jo('fitline_jo',X,Y,250);
ParMean = fitline_jo(X,Y);
ParStd  = std(Par);
