function [P1,P2]=anovacomp(x,c)
% ANOVACOMP(x,c) -- Performs an ANOVA on x (with each group having its own column),
%                   and a complex comparison using the weights listed in c
%
% NOTE: use 'NAN' to fill in non-values in the matrix x, DO NOT USE ZEROS!

% get number of columns
nCol = size(x,2);

% get number of rows
nRow = size(x,1);

% get size of c
sc = size(c);

% confirm c is a single-dimensioned array and of the right length
if sc(1) ~= 1 & sc(2) ~= 1
    error(' --> coefficients must be represented as a single-dimensioned array/list')
end
c = c(:)';
if length(c) ~= nCol
    error(' --> number of coefficients must equal the number of columns in data')
end

% get lengths of each column, disregrading NAN entries
lCol = sum(~isnan(x));

% get column means and psi
for cnt = 1:nCol
    colMean(cnt) = mean(x(find(~isnan(x(:,cnt))),cnt)); % we cannot calculate a mean if there is a NAN value in the list, thus we remove them first
end
psi = sum(colMean.*c);

% do standard ANOVA and get MSw and degrees of freedom
[P1, table, stats] = anova1(x,[],'off');
MSw = table{3,4};
df1 = table{2,3};
df2 = table{3,3};

% get F-statistic for complex comparison
F = psi.^2 / ...
    ( MSw .* sum( c.^2 ./ lCol ) );

% get P value
P2 = 1 - fcdf(F,df1,df2);

disp(' ')
disp([' Standard ANOVA P-value            = ' num2str(P1,'%06.4f')]);
disp([' Complex Comparison ANOVA P-value  = ' num2str(P2,'%06.4f')]);
disp(' ')

return