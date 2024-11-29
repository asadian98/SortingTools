function imagesc_plot(Ttime,datain,plottitle,scale,marks1,marks2);
% Generic plotting function for any analog data.
% Call: imagesc_plot(Ttime,datain,plottitle,scale,marks1,marks2)
%   Ttime: time scale. Must match the columns of Data
%   Data: each row is a different trial
%   title: title of subplot
%   scale: minTime, maxTime, maxEMG, minEMG (4th variable is optional; 0 if not specified)
%   marks1: optional set of marks, relative to 0 of Ttime. NaNs not plotted
%   marks2: optional set of marks, relative to 0 of Ttime. NaNs not plotted

if ~isempty(datain)
    if size(datain,1) > 1
        
        minTime = scale(1);
        maxTime = scale(2);
        maxval = scale(3);
        
        if length(scale) == 4
            minval = scale(4);
        else
            minval = 0;
        end
        
        
        hold on;
        title(plottitle)
        imagesc(Ttime,1:size(datain,1),datain,[minval maxval]);
        colorbar
        axis([minTime maxTime 1 size(datain,1)]);
        line([0 0],[0 size(datain,1)]);
        if nargin == 6
            for i = 1:length(marks1) 
                
                if(~isnan(marks1(i)))
                %        plot(marks1(i),i,'ws')
                plot(marks1(i),i,'ko','MarkerFaceColor','w')
                end
            end
            for i = 1:length(marks2)
                if(~isnan(marks2(i)))
                plot(marks2(i),i,'wo','MarkerFaceColor','r')
                end
            end
        elseif nargin == 5
            for i = 1:length(marks1)
                if(~isnan(marks1(i)))
                %        plot(marks1(i),i,'ws')
                plot(marks1(i),i,'ko','MarkerFaceColor','w')
                end
            end
        end
    end
end
