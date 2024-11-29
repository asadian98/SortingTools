% [Par]=fitline_jo(x,y)
%    fit a straight line through data set (x,y)
%
%    Par=[Gain,Offset]
%    
%    Jeroen Goossens

function par=fitline_jo(x,y)

par=polyfit(x,y,1);