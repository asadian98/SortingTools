function ROCarea = calcROC_BC(A,B,plotflag)
% Function spits out ROC area for distribution, stepping through 100 intermediate values between minimum and maximum
% values of A and B
%tic

C = A;
A = B;
B = C;


minint = floor(min(min(A),min(B))); maxint = ceil(max(max(A),max(B)));
stepsize = (maxint-minint)/500;
if nargin <= 2; plotflag = 0; end


ROCplot = [];
for t = minint:stepsize:maxint
    ROCplot = [ROCplot; sum(B>=t)/length(B) sum(A>=t)/length(A)];
end

% New code to fix bug
A = unique(ROCplot,'rows');
ROCarea = 0;
for i = 2:size(A,1)
	deltax = A(i,1) - A(i-1,1);
	deltay = (A(i,2) - A(i-1,2))/2 + A(i-1,2);
	ROCarea  = ROCarea + (deltax * deltay);
end
% 
% if ROCarea < 0.5
%     ROCarea = 1 - ROCarea;
% end

if plotflag == 1
    plot(ROCplot(:,1),ROCplot(:,2));
    text(0.6, 0.2,num2str(ROCarea))
end


%toc