%[p,d]=ks2test_jo(x1,y1,x2,y2)
%  Two dimensional Kolmogorov-Smirnov test. Tests wether
%  the two dimensional distributions (x1,y1) and (x2,y2)
%  are different. See also Press et.al. pp 645
%
%  p = significance level
%  d = K-S statistic
%
%  Jeroen Goossens

function [p,d]=ks2test_jo(x1,y1,x2,y2);

n1 = length(x1);
n2 = length(x2);

d1=0;
for j=1:1:n1
  % first use points in the first set as origins
  [fa,fb,fc,fd]=quadct_jo(x1(j),y1(j),x1,y1);
  [ga,gb,gc,gd]=quadct_jo(x1(j),y1(j),x2,y2);
  d1 = max([d1,abs(fa-ga)]);
  d1 = max([d1,abs(fb-gb)]);
  d1 = max([d1,abs(fc-gc)]);
  d1 = max([d1,abs(fd-gd)]);
end;

d2 = 0;
for j=1:1:n2
  % Then use points in the second set as origins
  [fa,fb,fc,fd]=quadct_jo(x2(j),y2(j),x1,y1);
  [ga,gb,gc,gd]=quadct_jo(x2(j),y2(j),x2,y2);
  d2 = max([d2,abs(fa-ga)]);
  d2 = max([d2,abs(fb-gb)]);
  d2 = max([d2,abs(fc-gc)]);
  d2 = max([d2,abs(fd-gd)]);
end;

% average the K-S statistics
d = 0.5*(d1+d2);   
sqen = sqrt( n1*n2/(n1+n2) );

% get linear correlation coefficient for each set
r1=pearsn_jo(x1,y1);
r2=pearsn_jo(x2,y2);
rr=sqrt(1-0.5*(r1*r1+r2*r2));

% estimate the probatility 
p = probks_jo( d*sqen/(1+rr*(0.25-0.75/sqen)) );


