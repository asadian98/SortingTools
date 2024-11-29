%[fa,fb,fc,fd]=quadct_jo(x,y,xx,yy)
%  Function used by ks2test
%
%  Jeroen Goossens

function [fa,fb,fc,fd]=quadct_jo(x,y,xx,yy)

%  Given an origin (x,y), and coordinates xx,yy this
%  procedure counts how many of the points are in each
%  quadrant around the origin, and returns the normalized
%  fractions. Quadrants are labeled alphabetically, 
%  counterclockwise from the upper right. 
%  See also Press et. al. pp 648

na = sum( (yy>y)  & (xx>x)  );
nb = sum( (yy>y)  & (xx>=x) );
nc = sum( (yy<=y) & (xx>=x) );
nd = sum( (yy<=y) & (xx>x)  );

nn = length(xx);
ff = 1/nn;
fa = ff*na;
fb = ff*nb;
fc = ff*nc;
fd = ff*nd;
