function outangle = rotaryangle(Ehamp, Evamp)
% function rotaryangle
% Designed by bdc on Feb 01/02 to return the rotary direction of vector [Ehamp Evamp], 
% starting from 90R and proceeding in a counterclockwise direction

if Ehamp == 0
    if Evamp > 0
        outangle = 90;
    elseif Evamp < 0
        outangle = 270;
    else
        outangle = -999;
    end
elseif Evamp == 0
    if Ehamp > 0
        outangle = 0;
    elseif Ehamp < 0
        outangle = 180;
    else
        outangle = -999;
    end
else
    outangle = atan(Evamp/Ehamp) * (360/(2*pi)); %Quadrant 1
    % Work four quadrants as necessary 
    if Ehamp < 0 & Evamp > 0    %Quadrant 2
        outangle = 180 + outangle;
    elseif Ehamp < 0 & Evamp < 0 % Quadrant 3
        outangle = 180 + outangle;
    elseif Ehamp > 0 & Evamp < 0 % Quadrant 4
        outangle = 360 + outangle;
    end
end