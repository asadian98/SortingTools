function onset_latency = find_onset_latency_ROC(time_re_target,ROC,visual_interval,criteria,plotflag)
% Function designed by BDC on Oct 7 2022 to find onset latency of mean
% spike density function for Puk Nuijten's analysis of visual onset latency
% when movements are only in one direction. Important point is that this
% function "detrends" the data between baseline_start and baseline_end,
% which is also the segment used to determine the mean plus the
%
% Outputs NaN if no value is found
%
% Inputs:
%   - time_re_target: the timeframe covered by the spike density function
%   - ROC: the time series ROC
%   - visual_interval: A 2 number array for where you are scanning for
%        potential onsets (e.g., [30 100]);
%   - criteria: A three number array, consisting of:
%        - the ROC value to exceed
%        - the critical number of points above threshold
%        - the number of points to scan for being above threshold
%       So [.6 8 10] means you find first point above a threshold of 0.6 where 8 of the next 10 points are above threshold
%   - plotflag. Set to 1 if you want to see plots during development. Defaulted to zero
%
% USAGE
%   onset_latency = find_onset_latency_ROC(ROC(:,1),ROC(:,2),[30 100],[.6 8 10],0)

threshold = criteria(1); numover = criteria(2); numover_scan = criteria(3);

onset_latency = NaN;    %Default output


if nargin < 5; plotflag = 0; end% By Default do not plot anything


if plotflag; h2 = figure; set(gcf,'Name','ROC Onset analysis');
    plot(time_re_target,ROC); hold on;
    plot([min(time_re_target) max(time_re_target)],[threshold threshold],'k--');
    axis([min(time_re_target) max(time_re_target) 0 1]);
end


flag_threshold = 0; 
start_point = find(time_re_target == visual_interval(1));end_point = find(time_re_target == visual_interval(2));
for i = start_point:end_point
    if ~flag_threshold && ROC(i) >= threshold
        if sum(ROC(i:i+numover_scan-1)>=threshold) >= numover
            onset_latency = time_re_target(i);
            flag_threshold = 1;
            if plotflag
            YAXIS = ylim; Ymax = YAXIS(2); Ymin = YAXIS(1); XAXIS = xlim; Xmin = XAXIS(1); Xmax = XAXIS(2);
            plot([onset_latency onset_latency],[0 Ymax],'m--')
            text(Xmin + 0.1*(Xmax-Xmin),Ymin + 0.1*(Ymax-Ymin), strcat('Onset latency = ',num2str(onset_latency)));
            end
        end
%     elseif(~flag_threshold && ROC(i) <= 0.5 + (0.5 - threshold))
%         if sum(ROC(i:i+numover_scan-1)<=0.5 + (0.5 - threshold)) >= numover
%             onset_latency = time_re_target(i);
%             flag_threshold = 1;
%             if plotflag
%                 YAXIS = ylim; Ymax = YAXIS(2); Ymin = YAXIS(1); XAXIS = xlim; Xmin = XAXIS(1); Xmax = XAXIS(2);
%                 plot([onset_latency onset_latency],[0 Ymax],'m--')
%                 text(Xmin + 0.1*(Xmax-Xmin),Ymin + 0.1*(Ymax-Ymin), strcat('Onset latency = ',num2str(onset_latency)));
%             end
%         end
    end
end

      
