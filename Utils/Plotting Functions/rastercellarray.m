function [outputcellarray] = rastercellarray(datain, i, NU, targetpos, trialtype, fid, nLFP_chans, LFPchan_names, adtodeg, event, timespan)
% function rastercellarray
% Designed by bdc on June 23/02 to produce a cell array for output for trials which
% match the Unit, target and trialtype combination. Outputcellarray is a 3 X n cell array, where n = number of trials
% Field 1 = spiketimes, zeroed to the pretime given in "timespan"
% Field 2 = Eh, zeroed to the event time and spanning the pre time given in "timespan"
% Field 3 = Ev, zeroed to the event time and spanning the pre time given in "timespan"
% datain is the entire datastr, i is the array of correct trials to scan

outputcellarray = [];
z = 0;  % Record counter for outputstructure

for b = 1:size(i,1) % Scan correct trials only
    a = i(b);
    if (datain(a).header(19) == trialtype) & (datain(a).header(8) == targetpos(1)) & (datain(a).header(9) == targetpos(2))
        z = z + 1;  % Increment record counter
        spiketimes = datain(a).tspikes{NU}; % Grabs spike times for desired unit in trial
        
        % Determine shift time of aligning event
        if event < 90
            shifttime = datain(a).events(event) + timespan(1);
        else
            sacnum = min(find(datain(a).saccades(:,4) > 80));   % Time of first non-anticipatory saccade after FP offset            
            shifttime = datain(a).saccades(sacnum,1) + timespan(1);
        end
        
        spiketimes = spiketimes - shifttime;   % Zeroes spiketimes around pretime relative to desired event
        outputcellarray{1,z} = spiketimes(spiketimes > 0 & spiketimes < (timespan(2) - timespan(1)));    % Grabs only spikes within window defined by timespan
                
        % Find eye traces during window defined by timespan
        pointer_pos = ((datain(a).tsegment(1) - 1) * nLFP_chans) * 2;
        fseek(fid, pointer_pos, 'bof');
        conttraces = fread(fid, [nLFP_chans datain(a).tsegment(3)], 'int16');

        % Convert ad units to degrees, using conversion factor found in configuration file. 
        %Filter with soft sgolayfilt since Dan's gain is cranked very high. Introduces a minimal phase shift
        for x = 1:nLFP_chans
            if sum(LFPchan_names(x,1:4) == 'AD17') == 4 | sum(LFPchan_names(x,1:4) == 'Heye') == 4; rawEh = conttraces(x,:) / adtodeg(1);   % corresponds to AD17
            elseif sum(LFPchan_names(x,1:4) == 'AD18') == 4 | sum(LFPchan_names(x,1:4) == 'Veye') == 4; rawEv = conttraces(x,:) / adtodeg(2);   % corresponds to AD18
            end
        end
        Eh = sgolayfilt(rawEh,3,11);
        Ev = sgolayfilt(rawEv,3,11);
        
        outputcellarray{2,z} = Eh(shifttime:shifttime+timespan(2)-timespan(1));
        outputcellarray{3,z} = Ev(shifttime:shifttime+timespan(2)-timespan(1));        
        
    end     
end
