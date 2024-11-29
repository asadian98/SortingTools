% q = outlier_jo(X,t);
%   returns an index vector q that is 0 for outliers in set X 
%   otherwise 1. t = threshold, usually 3 to 4;
%   uses standardized value criterion ( see Neter et al, pp 82-83);

function q = outlier_jo(X,t);

if length(X) == 1;
    q = 1;
else
    % normalize data between -1 and 1
    mX = min(X);
    X  = X-mX;
    
    pX = max(X);
    X  = X/pX;
    X  = 2*X-1;
    
    % calculate outliers
    mX = mean(X);
    sX = std(X);
    
    N  = abs((X-mX)/sX);
    q = N < t;
end
