function [sacarray] = sacdetector(dEh,dEv,Eh,Ev,thresh,DiodeTime)
critical_thresh = thresh(1); start_thresh = thresh(2); end_thresh = thresh(3);

absvel=sqrt((dEh.^2) + (dEv.^2));
i = size(dEh, 2);

C = find(absvel > critical_thresh); %Indices above critical threshold 

sacarray = [];
if isempty(C);sacarray = []; trajectories = [];  %Detect saccades etc only if there is at least one point above the threshold...
else
    D = C(diff(C) > 1); %Yields indices above threshold where the differences between adjacent points is greater than one. One such point per saccade
    D(1,size(D,2) + 1) = max(C); %Gets value for last saccade in trial
    b = 0;
    
    for y = 1:size(D,2)
        if D(y) > 10 & D(y) < i-10
            b = b+1;
            a = max(find(absvel(10:D(y)) < start_thresh));
            if isempty(a)
                sacarray(b,1) = 10;
            else
                sacarray(b,1) = 10 + a;
            end
            c = min(find(absvel(D(y):i-10) < end_thresh));
            if isempty(c)
                sacarray(b,2) = i -10;
            else
                sacarray(b,2) = c + D(y) + 2;
            end
        end
    end
    % Fill in entries to saccade table
    % First get relevant information from the events array for this trial
    Ton = DiodeTime;        % DUMMY VARIABLES HERE
    FPoff = DiodeTime;
    Vibon = DiodeTime;
    
    % Remove duplicates
    if(~isempty(sacarray))

        ZZ = sacarray(:,1);
        duplicaterows = find(diff(ZZ)==0);
        if ~isempty(duplicaterows);sacarray = removerows(sacarray,duplicaterows);end
        
        % For each entry in the saccade array, calculate saccade parameters
        B = zeros(size(sacarray,1),18);
        sacarray = [sacarray B];    %Preallocates sacarray matrix with zeros
        for b = 1:size(sacarray,1)
            sacarray(b,3) = sacarray(b,1) - Ton;    %SRT relative to Cue on
            sacarray(b,4) = sacarray(b,1) - FPoff;  % SRT relative to FP off
            sacarray(b,5) = sacarray(b,1) - Vibon;  % SRT relative to Vibon
            sacarray(b,6) = Eh(sacarray(b,1));    % Eh on
            sacarray(b,7) = Ev(sacarray(b,1));    % Ev on
            sacarray(b,8) = Eh(sacarray(b,2));    % Eh off
            sacarray(b,9) = Ev(sacarray(b,2));    % Ev off
            sacarray(b,10) = sacarray(b,8) - sacarray(b,6); % Eh amp
            sacarray(b,11) = sacarray(b,9) - sacarray(b,7); % Ev amp
            sacarray(b,12) = sqrt((sacarray(b,10) ^2) + (sacarray(b,11)^2)); %Vectorial amplitude
            sacarray(b,13) = sacarray(b,2) - sacarray(b,1); % Saccade duration
            sacarray(b,14) = rotaryangle(sacarray(b,10),sacarray(b,11)); % Saccade direction
            if dEh(1,sacarray(b,1)) > 0; sacarray(b,15) = max(dEh(1,sacarray(b,1):sacarray(b,2))); % Peak dEh
            else sacarray(b,15) = min(dEh(1,sacarray(b,1):sacarray(b,2))); end
            if dEv(1,sacarray(b,1)) > 0; sacarray(b,16) = max(dEv(1,sacarray(b,1):sacarray(b,2))); % Peak dEv
            else sacarray(b,16) = min(dEv(1,sacarray(b,1):sacarray(b,2))); end
            sacarray(b,17) = max(absvel(1,sacarray(b,1):sacarray(b,2))); % Peak vectorial velocity
            sacarray(b,18) = find(max(absvel(1,sacarray(b,1):sacarray(b,2)))); % Time to peak vel
            sacarray(b,19) = skewness(absvel(1,sacarray(b,1):sacarray(b,2))'); % Skewness of absolute vel profile
            
            
            % Save position and velocity trajectories
            % Access trajectory by: datastr(n).trajectories(b).Eh, where n = trial number, b = saccade number
            
            trajectories(b).Eh = Eh(1,sacarray(b,1):sacarray(b,2)); 
            trajectories(b).Ev = Ev(1,sacarray(b,1):sacarray(b,2)); 
            trajectories(b).dEh = dEh(1,sacarray(b,1):sacarray(b,2)); 
            trajectories(b).dEv = dEv(1,sacarray(b,1):sacarray(b,2)); 
            
            quarter = fix((sacarray(b,2)-sacarray(b,1))/4) + sacarray(b,1);
            qvectEh = Eh(1,quarter) - Eh(1,sacarray(b,1));
            qvectEv = Ev(1,quarter) - Ev(1,sacarray(b,1));
            sacarray(b,20) = rotaryangle(qvectEh, qvectEv);
            
        end

    end

end


return
