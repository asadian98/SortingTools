function [chsq,prob,df]=chsone(bins,ebins,knstrn)
%CHSONE	Chi-Square Significance Test
%	[CHSQ,PROB,DF]=chsone(BINS,EBINS,KNSTRN) performs a Chi-Square
%	Test on the data histogram in BINS to the distribution in
%	EBINS.  KNSTRN is the number of constraints used in the model
%	after the data was collected, usually zero.  CHSQ is Chi-Squared
%	value, and PROB is the "significance". A small value of PROB
%	indicates a significant difference between the distributions
%	BINS and EBINS.  DF is the number of degrees of freedom.
%
%	Translation of CHSONE from Numerical Recipes, 1st ed. p. 488. or
%	2nd ed. p. 621.
%	Michael Maurer, 14 May 1991

if nargin<3,
   knstrn=0;
end
if any(ebins<0),
   error('Non-positive expected number in ebins')
end
I=find((bins~=0) | (ebins~=0));		% omit 0=nj=Nj from sum
nbins=length(I);
df=nbins-1-knstrn;
chsq=sum(((bins(I)-ebins(I)).^2)./ebins(I));
%prob=1-gamma(df/2,chsq/2);
%prob=1-gammainc(df,chsq);
prob=1-chi2cdf(chsq,df);