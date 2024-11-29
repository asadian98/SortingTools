function [chance,observed] = discriminate_bc(A,B,iterations)
% Simple function to run discriminant analysis on A and B, repeated
% "iterations" number of times

if nargin < 3
    iterations = 1
end

RESULTS = [];

for i=1:iterations

    i
    
    RESULTS(i,1) = length(A)/( length(A) + length(B) );
    
    % Randomly shuffle input arrays
    AA = A(randperm(length(A))); BB = B(randperm(length(B))); 
    AA_mean = mean(AA(1:round(length(A)/2)));
    BB_mean = mean(BB(1:round(length(B)/2)));
    
    % put in a 1 for trial type A, 0 for trial type B
    Trial_array=[ones(1,length(round(length(A)/2)+1:length(A))) zeros(1,length(round(length(B)/2)+1:length(B)))];

    
    Classified_array = [];
    Chance_array = [];
   z = 0;
    
    Rand_indices = randperm(length(Trial_array));
    for x = round(length(A)/2)+1:length(A);
        % Classify based on whether osbservation from AA is closer to
        % mean(AA) or mean(BB). 
        z = z+1;
        if abs(AA(x) - AA_mean) < abs(AA(x) - BB_mean)  % If so, then observation is closer to A than B
            Classified_array(end+1) = 1;    % Classified as A
        else
            Classified_array(end+1) = 0;    % Classified as B
        end
        
        Chance_array(end+1) = Trial_array(Rand_indices(z));
        
        
    end

    for y = round(length(B)/2)+1:length(B);
        z = z+1;
        % Classify based on whether osbservation from BB is closer to
        % mean(AA) or mean(BB). 
        if abs(BB(y) - AA_mean) < abs(BB(y) - BB_mean)  % If so, then observation is closer to A than B
            Classified_array(end+1) = 1;         % Classifed as A
        else
            Classified_array(end+1) = 0;         % Classified as B
        end
        
        Chance_array(end+1) = Trial_array(Rand_indices(z));
    
    end

    Difference_array = abs(Trial_array - Classified_array);
    Difference_chance_array = abs(Trial_array - Chance_array);
    RESULTS(i,1) = mean(Difference_chance_array);
    RESULTS(i,2) = mean(Difference_array);
    
    
end

chance = mean(RESULTS(:,1);
observed = mean(RESULTS(:,2);


if nargin == 1
    display_opt = 0;
end

data_mean = mean(datain);
data_std = std(datain);
data_stderr = stderr_bc(datain);
data_n = length(datain);

if display_opt == 1;    % Echo text
    disp(strcat('Mean=',num2str(data_mean)));
    disp(strcat('Std=',num2str(data_std)));
    disp(strcat('Stderr=',num2str(data_stderr)));
    disp(strcat('n=',num2str(data_n)));
    
end
