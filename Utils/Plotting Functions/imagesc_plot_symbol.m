function imagesc_plot_symbol(Ttime,datain,plottitle,scale,marks1,marks2,symbols);
% Generic plotting function for any analog data.
% Call: imagesc_plot(Ttime,datain,plottitle,scale,marks1,marks2)
%   Ttime: time scale. Must match the columns of Data
%   Data: each row is a different trial
%   title: title of subplot
%   scale: minTime, maxTime, maxEMG, minEMG (4th variable is optional; 0 if not specified) 
%   marks1: optional set of marks, relative to 0 of Ttime
%   marks2: optional set of marks, relative to 0 of Ttime
%   symbols: color and types of symbols for marks 1 and marks

minTime = scale(1);
maxTime = scale(2);
maxval = scale(3);

if length(scale) == 4
    minval = scale(4);
else
    minval = 0;
end

if nargin < 7
    color1 = 'k';
    colorline1 = 'w';
    marker1 = 'o';
    color2 = 'w';
    colorline2 = 'k';
    marker2 = '>';
else
    color1 = symbols(1);
    colorline1 = symbols(2);
    marker1 = symbols(3);
    color2 = symbols(4);
    colorline2 = symbols(5);
    marker2 = symbols(6);
end


hold on;
title(plottitle)
imagesc(Ttime,1:size(datain,1),datain,[minval maxval]);
colorbar
axis([minTime maxTime 1 size(datain,1)]);
line([0 0],[0 size(datain,1)]);
if nargin == 6 | nargin == 7
    for i = 1:length(marks2)
        plot(marks2(i),i,strcat(colorline2,marker2),'MarkerFaceColor',color2)
    end
    for i = 1:length(marks1)
%        plot(marks1(i),i,'ws')
        plot(marks1(i),i,strcat(colorline1,marker1),'MarkerFaceColor',color1)
    end
elseif nargin == 5
    for i = 1:length(marks1)
%        plot(marks1(i),i,'ws')
        plot(marks1(i),i,strcat(colorline1,marker1),'MarkerFaceColor',color1)
    end
end
   