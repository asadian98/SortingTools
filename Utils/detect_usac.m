function [microsac_table] = detect_usac(Gh_all,Gv_all,dGh_all,dGv_all,Ttime)
% Function designed by BDC on Feb 13 2014 to try to detect microsaccades
% during the provided timeframe.
% Works in a fairly simple fashion to use a fairly low velocity criteria of
% ~25deg/s and crossings of 10deg/s.
% NOTE THIS IS NOT MEANT TO REPLACE INSPECTION -- THIS IS MAINLY USED AS A
% FIRST-PASS AUTODETECTOR FOR 2014 CIHR GRANT

% CONSTRUCTION OF MICROSAC_TABLE.
% Each row is a different usac
% Col 1: trial number
% Col 2: onset time in trial re Target, based on Ttime
% Col 3: offset time in trial, based on Ttime
% Col 4: trial type indicator (-1 = cue left, 0 = no cue, 1 = cue right)
% Col 5: Ghamp
% Col 6: Gvamp
% Col 7: Vectorial G amp
% Col 8: Direction, 0 = straight right, proceeding counterclockwise
% Col 9: dGh peak
% Col 10: dGv peak
% Col 11: dG_peak
% Col 12: duration

microsac_table = [];

thresh = [25 10 10];    % detected sac must exceed thresh(1), onset and offset marked at thresh(2) and thresh(3)
Gcritical_thresh = thresh(1);
Gstart_thresh = thresh(2);
Gend_thresh = thresh(3);


for x = 1:size(Gh_all,1)
%     x
    dGh = dGh_all(x,:); dGv = dGv_all(x,:); Gh = Gh_all(x,:); Gv = Gv_all(x,:);
%     figure; subplot(2,1,1); hold on;plot(Ttime,Gh,'r-');plot(Ttime,Gv,'g-'); subplot(2,1,2); hold on; plot(Ttime,dGh,'r-');plot(Ttime,dGv,'g-');plot(Ttime,sqrt((dGh.*dGh) + (dGv.*dGv) ) ,'k.')
    
%    absvel=sqrt((dGh.^2) + (dGv.^2));
    absvel=abs(dGh);
    i = size(dGh, 2);
    
    C = find(absvel > Gcritical_thresh); %Indices above critical threshold
    
    
    if ~ isempty(C)%Detect movements etc only if there is at least one point above the threshold...
        D = C(diff(C) > 25); %Yields indices above threshold where the differences between adjacent points is greater than twenty-five. One such point per movement
        D(1,size(D,2) + 1) = max(C); %Gets value for last movement in trial
        b = 0;
        cursize = size(microsac_table,1);
        
        for y = 1:size(D,2)
            if D(y) > 10 & D(y) < i-10
                b = b+1;
                a = max(find(absvel(10:D(y)) < Gstart_thresh));
                if isempty(a)
                    microsac_table(end+1,1) = x;
                    microsac_table(end,2) = min(Ttime)+10;
                else
                    microsac_table(end+1,1) = x;
                    microsac_table(end,2) = min(Ttime)+10+a;
                end
                c = min(find(absvel(D(y):i-10) < Gend_thresh));
                if isempty(c)
                    microsac_table(end,3) = i-10+min(Ttime);
                else
                    microsac_table(end,3) = c + D(y) + 2+min(Ttime);
                end
            end
        end
        % For each entry in the movement array, calculate movement parameters
        % Remove duplicates
        for z = 1:b
            onset = find(Ttime == microsac_table(cursize+z,2));
            offset = find(Ttime == microsac_table(cursize+z,3));
            microsac_table(cursize+z,4) = 0;    % Indicating FP
            microsac_table(cursize+z,5) = Gh(offset)-Gh(onset);    % Ghamp
            microsac_table(cursize+z,6) = Gv(offset)-Gv(onset);    % Gvamp
            Ghamp = Gh(offset)-Gh(onset);    % Ghamp
            Gvamp = Gv(offset)-Gv(onset);    % Gvamp
            microsac_table(cursize+z,7) = sqrt( (Ghamp*Ghamp)+(Gvamp*Gvamp));    % Vectorial amp
            microsac_table(cursize+z,8) = rotaryangle(Ghamp,Gvamp);
            if Ghamp >=0; microsac_table(cursize+z,9) = max(dGh(onset:offset)); else; microsac_table(cursize+z,9) = min(dGh(onset:offset));end
            if Gvamp >=0; microsac_table(cursize+z,10) = max(dGv(onset:offset)); else; microsac_table(cursize+z,10) = min(dGv(onset:offset));end
            microsac_table(cursize+z,11) = max(absvel(onset:offset));
            microsac_table(cursize+z,12) = offset-onset;
        end
    end
    
    
end

return

