function ttest

%program to determine difference of means for Stats assign 1
%created by BC on Jan21



%Initialize arrays
linux = [691 535 634 733 672 668 684 518 505 697 677 623 539 551 554 667 673 519 601 659 664 673 511 602 654];
MS = [653 505 482 617 536 612 621 543 464 580 559 497 592 515 547 494 483 522 452 571 492 483 526 455 570];

% Histograms and stats
maxRT = ceil(max(max(linux),max(MS))/10) * 10;

%Within samples t test
Di = (linux)' - (MS)'; 
D1 = sum(Di);
Dii = Di .* Di;
D2 = sum(Dii);
n = numel(linux);
DF = n - 1;
t = D1/sqrt(((n*D2)-(D1 * D1))/ (n-1))



%standard error of the mean
SEms = std(MS)/sqrt(n);
SElinux = std(linux)/sqrt(n);




EDGES = 450:10:maxRT;
N_linux = histc(linux,EDGES)/length(linux);
N_MS = histc(MS,EDGES)/length(MS);





figure;
subplot (2,2,1);hold on;
title ('Linux latencies Independant samples')
xlabel('Secs')
ylabel('%')
bar(EDGES,N_linux,'histc','w');
text(425,0.15,strcat('Mean = ',num2str(mean(linux)),'+/-',num2str(std(linux)),'. n=',num2str(length(linux))))
[H,significance,ci,stats]=ttest2(linux,MS);
text(425,0.13,strcat('t = ',num2str(stats.tstat), '  DF = ',num2str(stats.df)))
text(425,0.11,strcat('P (Linux vs MS) =',num2str(significance))) 
text(425,0.09, strcat('95% C.I. = ',num2str(ci)))
axis([400 750 0 0.17])

subplot (2,2,2);hold on;
title ('MS latencies Independant Samples')
xlabel('Secs')
ylabel('%')
bar(EDGES,N_MS,'histc','w');
text(550,0.14,strcat('Mean = ',num2str(mean(MS)),'+/-',num2str(std(MS)),'. n=',num2str(length(MS))))
axis([400 750 0 0.17])

subplot (2,2,3);hold on;
title ('Linux latencies Within Samples')
xlabel('Secs')
ylabel('%')
bar(EDGES,N_linux,'histc','w');
text(425,0.15,strcat('Mean = ',num2str(mean(linux)),'+/-',num2str(std(linux)),'. n=',num2str(length(linux))))
text(425, 0.13,strcat('t = ', num2str(t), '   Alpha level = 0.01  DF = ',num2str(DF) ))
text(425, 0.11,strcat('95% C.I. =  55.767     114.873'))
axis([400 750 0 0.17])



subplot (2,2,4);hold on;
title ('MS latencies Within Samples')
xlabel('Secs')
ylabel('%')
bar(EDGES,N_linux,'histc','w');
text(450,0.14,strcat('Mean = ',num2str(mean(MS)),'+/-',num2str(std(MS)),'. n=',num2str(length(MS))))
axis([400 750 0 0.17])



return