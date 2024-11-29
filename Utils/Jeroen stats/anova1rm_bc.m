function [F,P] = anova1rm_bc(table)
% Function created by bdc on July 15th to perform one-way RM ANOVA on table
% Table is arranged so that the rows are the subjects, and the columns
% are the treatments.
% Algorithm based on Zar, 1999, Biostatistical Analysis, 4ed, pg 255-259
% Post-hoc Tukey test performed by building STATS structure as required for
% MULTCOMPARE function. NOT SURE THIS IS KOSHER!


SUMXij = sum(sum(table));
table2 = table .* table;
SUMX2ij = sum(sum(table2));

C = (SUMXij * SUMXij)/(size(table,1)*size(table,2));

TotalSS = SUMX2ij - C;

rowsum = sum(table,2);
sumrowsum2 = sum(rowsum .* rowsum);
SubjectsSS = sumrowsum2/size(table,2) - C;

WithinSubjectsSS = TotalSS - SubjectsSS;

colsum = sum(table,1);
sumcolsum2 = sum(colsum .* colsum);
TreatmentSS = sumcolsum2/size(table,1) - C;

RemainderSS = WithinSubjectsSS - TreatmentSS;

TotalDF = size(table,1) * size(table,2) - 1;
SubjectsDF = size(table,1) - 1;
WithinSubjectsDF = TotalDF - SubjectsDF;
TreatmentDF = size(table,2) - 1;
RemainderDF = WithinSubjectsDF - TreatmentDF;

TreatmentMS = TreatmentSS / TreatmentDF;
RemainderMS = RemainderSS / RemainderDF;

F = TreatmentMS / RemainderMS;
P = 1 - fcdf(F,TreatmentDF,RemainderDF);

lineA = 'Analysis of Variance Table';
lineA2 = '';
lineB = sprintf('%s\t%s\t\t\t%s\t\t%s\t','Source of var','SS','DF','MS');
lineC = sprintf('%s\t\t\t%.2f\t%d\t','Total',TotalSS,TotalDF);
lineD = sprintf('%s\t\t%.2f\t%d\t','Subjects',SubjectsSS,SubjectsDF);
lineE = sprintf('%s\t\t%.2f\t\t%d\t\t%.2f\t','Treatments',TreatmentSS,TreatmentDF,TreatmentMS);
lineF = sprintf('%s\t\t%.2f\t\t%d\t\t%.2f\t','Remainder',RemainderSS,RemainderDF,RemainderMS);
lineG = '';
lineH = sprintf('%s\t%.2f\t',strcat('F(1,',num2str(TreatmentDF),',',num2str(RemainderDF),') = '),F);
lineI = sprintf('%s\t%.6f\t','P = ',P);

sprintf('%s\n',lineA,lineA2,lineB,lineC,lineD,lineE,lineF,lineG,lineH,lineI)

% Create structure for multiple comparisons test

for i = 1:size(table,2)
    STATS.gnames(i,1) = num2str(i);
    STATS.n(i) = size(table,1);    
    STATS.means(i) = mean(table(:,i));
end
STATS.source = 'anova1';
STATS.df = RemainderDF;
STATS.s = RemainderMS;

COMPARISON = multcompare(STATS,0.05,'on','tukey-kramer')

return

