function p_corr = holm_bonferroni(pvals, alpha)
% Finds the holm-bonferroni corrected pval for a group of unsorted pvals.
% See wikipedia page on this topic.
%
% Algorithm first sorts all Pvals from lowest to highest, and for each
% p-val tests whether it is lower than the Bonferroni corrected value for
% the remaining number of tests. Keeps repeating this until it finds a test
% that fails. This then is the returned corrected p-value.
%
% By default, alpha = 0.05
if nargin < 2
    alpha = 0.05;
end

if size(pvals,1) < size(pvals,2)
    pvals = pvals';
end
n = length(pvals);
p_sorted = sort(pvals,1,'descend');
p_sorted = cat(2,p_sorted,alpha./[1:n]');

crit = min(find(p_sorted(:,1)<p_sorted(:,2)));
if isempty(crit) | isnan(crit)
    p_corr = alpha/n; %Worst case scenario is straight Bonferroni corrected val
else
    p_corr = alpha/crit;
end

return





