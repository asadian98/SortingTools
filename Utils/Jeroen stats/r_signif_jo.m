% p = r_signif_jo(R,N)
%   calculte two-side significance level of a 
%   correlation coefficient R based on N data 
%   points. (N>=10)
% 
%   Jeroen Goossens

function p = r_signif_jo(R,N)

% significance level of a linear correlation
% see: Press et. al.
%      Numerical Recipes 2ed edition p 636

df = N-2;                              % degrees of freedom 
t  = R .* sqrt( df./(1-R.*R) );       % t-value 

x  = df./(df+t.*t);
a  = df./2;
b  = 1/2;   
if x > 1 & x < 1.29
    x
    disp('Adjusting x val to 1')
    p = betainc(0.9999,a,b); % ALLOWING SOME WIGGLE ROOM?
else
    p  = betainc(x,a,b);                       % p-value
end

% 
% function p_F = r_signif(R,N);

% in_par_len = 2;
% R_sq = R*R;
% F = [(in_par_len - 1), (N - 1 - in_par_len), ...
%     ((R_sq/(in_par_len - 1))/ ((1 - R_sq)/(N - 1 - in_par_len)))];
% p_F = beta(F(2)/(F(2) + F(1)*F(3)), F(2)/2, F(1)/2)
