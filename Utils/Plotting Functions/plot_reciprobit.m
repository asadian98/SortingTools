function plot_reciprobit(RT,plotcolor,plot_name);
% simple plotting function to plot RT data in reciprobit fashion
if nargin < 3; plot_name = 'Reciprobit plot';end
if nargin < 2; plotcolor = 'k';end


% SHADLEN'S cumprob function
RT_sorted = sort(RT);
p = cumsum(RT./RT)/(length(RT)+1);

rtmin = 10; rtmax = 1000; n = 2000;
rtax = [floor(.5*min(RT)):ceil(2*max(RT))];

plot(-1./RT_sorted, norminv(p), 'ko','MarkerSize',4,'MarkerFaceColor',plotcolor);
hold on;

for i = floor(min(RT)/100)*100:100:ceil(max(RT)/100)*100
    plot([-1/i -1/i],[-2.5 3],'b--')
    text(-1/i,-2.75,num2str(i))
end

return

