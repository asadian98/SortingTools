function [rasters] = raster(datain, subdir, convfnct, event, timespan)
% function raster
% Syntax: function [rasters] = raster(datain, subdir, convfnct, event, timespan)
% e.g., : [rasters] = raster('Dan0216b','Dan_vib',{'gaussian', 'sigma'}, 4, [-100 500])
% Designed by bdc on June 21/02 to print out rasters, average spike density and eye movements
% for all target positions, aligned on the desired event.
%   convfnct = Cell array of convolving function type (cell 1) and parameters (cell 2..n)
%   event = Aligning event (see events matrix). 4 = Cue on, 11 = FP off, 99 = First sac after FP off
%   timespan = array of [pretime and posttime], defining range to display around aligning event
%
% STRUCTURE OF RASTERS MULTI-DIMENSIONAL CELL ARRAY
%
% Dim#  Description
% 1     signal number (e.g., 'sig001a')
% 2     Target number (e.g., [-7.07 7.07]). Order set by "unique" function
% 3     Trial condition number (e.g., 90). Order set by "unique function
% 4     DATA. The data is organized in two cells
%  {1}   A 3xn cell array, where 1 = spike events, 2 = Eh and 3 = Ev for "n" trials, all zeroed to desired event
%  {2}   A 3 x timespan array containing the averaged spike waveform (1st row), and 2 X stddev (2nd and 3rd row)

% Load data and binary file

% Step 1, access and load data file and binary file
filename = strcat('d:\data\', subdir, '\', datain,'.mat');
load(filename);
fid = fopen(strcat('d:\data\', subdir, '\', datain, '.dat'),'r'); % See analyze_trial function for how binary data is accessed

% Step 2. Create targarray of target locations, and initialize rasters cell array

numunits = size(tspikes,2)
for i = 1:numunits
    unitarray(i,1:64) = tspikes{1,i};
end
unitarray = char(unitarray(:,1:7)); % Names of recorded signal

targarray = unique(headers(:,8:9),'rows');
numtargs = size(targarray,1)
xvals = unique(targarray(:,1),'rows'); xrange = (xvals(length(xvals)) - xvals(1,1)) * 1.15; xadd = min(xvals); % Values for plotting rasters (adds 15%)
yvals = unique(targarray(:,2),'rows'); yrange = (yvals(length(yvals)) - yvals(1,1)) * 1.15; yadd = min(yvals);
eyemax = max(max(xvals),max(yvals));eyemin = min(min(xvals),min(yvals));

condarray = unique(headers(:,19),'rows');
numcond = size(condarray,1)

rasters = cell(numunits,numtargs,numcond,2);  

% Default steps for convolving function
if convfnct{1} == 'gaussian'
   D = gausfunct(convfnct{2}); % produces row vector gaussian         
   gauswidth = convfnct{2};
end

% Compute other convovling functions here if required

% Step 3. Assign trial data and averaged data to 4th dimension

% Grab only successful trials, and get data for appropriate combination of unit, target, trial
i = find(headers(:,21));

for NU = 1:numunits
    maxspden(NU) = 0;
    for NTG = 1:numtargs
        for NCOND = 1:numcond 

            H = [];B = []; C = [];
            
            disp(strcat('Assigning data for: ', unitarray(NU,:), ' at position: ', num2str(targarray(NTG,:)), ' and trialtype:', num2str(condarray(NCOND))))
            
            rasters{NU,NTG,NCOND,1} = rastercellarray(datastr,i,NU,targarray(NTG,:),condarray(NCOND), fid, nLFP_chans, LFPchan_names, adtodeg, event, timespan);
                        
            % Create average spikedensity function here
            % First produce spike density function for each trial
            for TRNUM = 1:size(rasters{NU,NTG,NCOND,1},2)     
                nspikes = size(rasters{NU,NTG,NCOND,1}{1,TRNUM},2);   % Number of spikes in current  

                B = ceil(rasters{NU,NTG,NCOND,1}{1,TRNUM}); % Spike times rounded to nearest ms

                C = zeros(1,timespan(2) - timespan(1) + 1);   %Row of zeros for timespan under consideration
                C(B) = 1;   %Replaces 0 with 1 at spiketimes
               
                if sum(C) > 0    % If there are spikes, then convolve with function
                    if convolve.type == 'gaussian'
                        E = conv(C',D(:,2)') * 1000;   %Convolves spiketimes with gaussian, converting to spikes/s
                        H = [H; E(3*gauswidth : length(E) - 3*gauswidth - 1)'];   %Shift spden function back 3*gauswidth, and leave out last 3*gauswidth points (recall gaussian was defined for 6 * sigma)                          
                    end                       
                    % ENTER OTHER CONVOLVING FORMULAS HERE                       
                else % If there are no spikes, then spike density for current trial equals zero
                    H = [H; zeros(1,timespan(2) - timespan(1) + 1)];
                end                                    
            end
            I = mean(H); J = 2 * std(H);            
            rasters{NU,NTG,NCOND,2} = [I; I+J; I-J];    % 3 x timespan array with mean (row 1) and 2 stddev traces (rows 2 and 3)            
            if max(I) > maxspden(NU); maxspden(NU) = max(I);end %Pulls out maximum spike density value for this unit, for plotting.
        end
    end
end

% Step 5: plot data
Colortype = ['b', 'r', 'g', 'm', 'k'];   % Color types for different trial types, will cycle through if there are more than 5

for NU = 1:numunits

    maxyval = 20 * ceil(maxspden(NU)/20)
    
    figure; %new plot for each isolated unit
    set(gcf,'Position',[50 50 1180 900]); % Values in screen pixels

%    hold on;
    
    for NTG = 1:numtargs
        for NCOND = 1:numcond 

            xvalue = targarray(NTG,1); yvalue = targarray(NTG,2); xbase = (xvalue - xadd) / xrange + 0.03; ybase = (yvalue - yadd) / yrange + 0.03;

            numtrials = size(rasters{NU,NTG,NCOND,1},2);     % Number of successful trials of unit/target/condition combination
            increment = maxyval / numtrials;    % Base increment depends on number of trials and maximum value for this unit

            
           subplot('Position', [xbase (ybase+0.06) 0.08 0.03 ]); % Plot rasters
           % Extract raster values
            
            svals =[]; ypos = [];
                        
            for TRNUM = 1:numtrials     % For all trials
                
                yoffset = (numtrials * increment) - (TRNUM * increment);
                svals = [svals rasters{NU,NTG,NCOND,1}{1,TRNUM}];
                ypos = [ypos ones(1,size(rasters{NU,NTG,NCOND,1}{1,TRNUM},2)) * yoffset];
                      
            end

            h1 = plot(svals,ypos, 'Color',Colortype(NU),'MarkerSize',2,'Marker','o','LineStyle','none');
            axis ([0 (timespan(2) - timespan(1)) 0 maxyval]);                        
            set(gca,'XTick',0:200:(timespan(2) - timespan(1)),'XTickLabel',{(timespan(1):200:timespan(2))},'Fontsize',8);            


            
            subplot('Position', [xbase (ybase+0.03) 0.08 0.03 ]); % Plot spike density function
            h2 = plot(0:(timespan(2)-timespan(1)),rasters{NU,NTG,NCOND,2}(1,:),'Color',Colortype(NU),'LineWidth',3);
            axis ([0 (timespan(2) - timespan(1)) 0 maxyval]);                                    
            set(gca,'XTick',0:200:(timespan(2) - timespan(1)),'XTickLabel',{(timespan(1):200:timespan(2))},'Fontsize',8);

            
            subplot('Position', [xbase ybase 0.08 0.03 ]);  % Plot Eye movements
            hold on;
            for TRNUM = 1:numtrials     % For all trials
            
                plot(rasters{NU,NTG,NCOND,1}{2,TRNUM},'-b'); % Plots Eh
                plot(rasters{NU,NTG,NCOND,1}{3,TRNUM},'-r'); % Plots Ev
                
            end
            axis ([0 (timespan(2) - timespan(1)) eyemin eyemax]);                                    

            
        end
    end
end

%posarray


% APPENDED TEXT FROM SHMAPLAT

% Try to set colorbar at specific user-defined values


%subplot('Position',[0.15 0.62 0.25 0.2]);mapsclat('./data/OCIonlat','./data/OCIonND');
%caxis('manual');caxis([-28 -8]);colorbar
