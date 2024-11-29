function Y = pointfilter(X,halfpoints)
% Simple filter program designed by bdc to average each point in X over the
% desired number of points to produce Y
% Filters the rows, i.e., each row is a different trial

window = halfpoints * 2;

if halfpoints == 0;
    Y = X;
    return
end

Y = zeros(size(X,1),size(X,2));

for i = 1:size(X,1)
    Z = X(i,:);
    for j = 1:halfpoints
        Y(i,j) = mean(Z(j:j+window));
%        Y(i,j) = (sum(Z(j:j+window))) / (window+1);
    end
    
    for k = size(X,2)-halfpoints:size(X,2)
        Y(i,k) = mean(Z(k-window:k));
%        Y(i,k) = (sum(Z(k-window:k))) / (window+1);
    end

    for j = halfpoints+1:size(X,2)-halfpoints-1
        Y(i,j) = mean(Z(j-halfpoints:j+halfpoints));
%        Y(i,j) = (sum(Z(j-halfpoints:j+halfpoints))) / ((halfpoints*2)+1);
    end
end
