%[p,d]=kstest_jo(data1,data2)
%  Kolmogorov-Smirnov test. Test whether the data1 
%  and data2 are drawn from different distributions
%  Small values of p indicate significant difference
%  See Press et. al. pp 623
%  
%  p = signifficance level
%  d = K-S statistic 
%
%  Jeroen Goossens

function [p,d]=kstest_jo(data1,data2)

N1 = length(data1);
N2 = length(data2);

data1 = sort(data1);
data2 = sort(data2);

j1  = 1;  j2  = 1;
fn1 = 0;  fn2 = 0;
d   = 0;

while (j1<=N1) & (j2<=N2),  % if we are not done
  d1 = data1(j1);
  d2 = data2(j2);
  if (d1<=d2),              % next step is in data1
    fn1 = j1/N1;
    j1  = j1+1;
  end;
  if (d2<=d1),              % next step is in data2 
    fn2 = j2/N2;
    j2  = j2+1;
  end;
  dt = abs(fn2-fn1);
  if (dt>d),
    d = dt;
  end;
end;

% probability  is good enough if N>=4
Ne = (N1*N2)/(N1+N2);
p  = probks_jo( (sqrt(Ne) + 0.12 + 0.11/sqrt(Ne) ) * d );
