function dataout = fliparray(datain)
% Simple function to generate flipped array
% ie [1 2 3] becomes [3 2 1]

for i = 1:length(datain)
    dataout(i) = datain(length(datain)-i+1);
end