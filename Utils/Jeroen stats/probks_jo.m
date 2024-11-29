%p=probks_jo(alam)
%  Kolmogorov-Smirnov probability function
%  See Press et. al. pp 493
%
%  Jeroen Goossens
  
function p=probks_jo(lambda)

j=1:1000;

p=2*sum( (-1).^(j-1).*exp(-2*j.^2.*lambda.^2) );

if lambda<1/100, p=1; end;




