function ROCarea = calcROC(A, B)
% Function spits out ROC area for distribution, stepping through 100 intermediate values between minimum and maximum
% values of A and B
% Assumes that the hypothesis being tested is that A will be greater than B
% (ie, ROC area will be greater than 0.5 for divergent distributions).
%tic
A = sort(A);
B = sort(B);
if isequal(A,B);ROCarea = 0.5; return; end
minint = floor(min(A(1),B(1)));
maxint = ceil(max(A(end),B(end)));
stepsize = (maxint-minint)/100;

ROCplot = [];
for t = minint:stepsize:maxint
    ROCplot = [ROCplot; sum(B>=t)/length(B) sum(A>=t)/length(A)];
end
plot(ROCplot(:,1),ROCplot(:,2),'s-');
axis([0 1 0 1])

% New code to fix bug
A = unique(ROCplot,'rows');
ROCarea = 0;
for i = 1:size(A,1);
    if i == 1; 
        deltax = A(i,1); 
        deltay = A(i,2)/2; 
    else
        deltax = A(i,1) - A(i-1,1);
        deltay = (A(i,2) - A(i-1,2))/2 + A(i-1,2);
    end
    ROCarea  = ROCarea + (deltax * deltay);
end

%ROCplot(:,3) = [diff(ROCplot(:,2)) * -1; 0];
%ROCarea = sum(ROCplot(:,1) .* ROCplot(:,3));

%toc