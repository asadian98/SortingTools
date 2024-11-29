% [F,p,df]=sfanova_jo(X)
%  single factor anova for unequal sample sizes. 
%  Empty entries in X should be denoted as NaN.
%  Each row of X is an observation, each colum 
%  of X is a variable
%
%  Jeroen Goossens
  
function [F,p,df]=sfanova_jo(X)

Nrow = size(X,1);
Ncol = size(X,2);

% grand mean
Xgm  = mean(X(isnan(X)==0));

% means of columns
for i=1:Ncol,
  Xi(i) = mean( X(isnan(X(:,i))==0,i) );
  Si(i) = std ( X(isnan(X(:,i))==0,i) );
end;

% Variation between columns  due to differences 
% between column means Xi
Ni     = sum(isnan(X)==0);
SScols = sum(Ni.*((Xi-Xgm).^2));
DFcols = Ncol-1;
MScols = SScols/DFcols;

% Resedual variation due to differences between
% observations Xit and column means Xi
for i=1:Ncol,
  SS(i) = sum( (X(isnan(X(:,i))==0,i)-Xi(i)).^2 );
end;
SSres = sum(SS);
DFres = sum(Ni-1);
MSres = SSres/DFres;

SStot = sum( (X(isnan(X)==0)-Xgm).^2 );
DFtot = sum(Ni)-1;

% F-ratio
F  = MScols/MSres;
df = [DFcols,DFres];

% significance level
Xb = DFres/(DFres + DFcols*F);  % see press et al, pp 228 eqn 6.4.9
Ab = DFres/2; 
Bb = DFcols/2 ;
p  = 2*betainc(Xb,Ab,Bb);
if (p>1), p = 2-p; end;

fid=1;
fprintf(fid,'\n Sample Data \n');
fprintf(fid,'       n      mean       std\n');
for i=1:Ncol,
  fprintf(fid,'%8d  %8.2f  %8.2f\n',Ni(i),Xi(i),Si(i));
end;
fprintf(fid,'\n Analysis of Variance \n');
fprintf(fid,' Source       df        SS        MS         F         p\n');
fprintf(fid,' Factor %8d  %8.2f  %8.2f  %8.2f  %8.5f\n',DFcols,SScols,MScols,F,p);
fprintf(fid,' Error  %8d  %8.2f  %8.2f\n',DFres,SSres,MSres);
fprintf(fid,' Total  %8d  %8.2f\n',DFtot,SStot);