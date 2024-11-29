% R=linfit_jo(x,y)
%   linear regression y =ax+b
%   the result R = [ a std(a) b std(b) r p n]  

function R=linfit_jo(x,y);

r   = corrcoef(x,y);
r   = r(1,2);
n   = length(x);
p   = r_signif_jo(r,n);
[par stdv] = linefit_jo(x,y);
R  = [ par(1) stdv(1) par(2) stdv(2) r p n];
