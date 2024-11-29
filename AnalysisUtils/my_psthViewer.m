function my_psthViewer(Z_Spikes, Z_Analog, eventTimes, trGroups, window, plotclassification, plotdelayed, plotmap, varargin)
% function psthViewer(spikeTimes, clu, eventTimes, window, trGroups[, params])
%
% Controls: see the command window text on startup.
% params can include 'groupingName' (a string) and 'groupingLegend' (a cell
% array of strings)
%
% TODO:
% XX if plotting all traces, color raster and sort
% - don't replot traces just change X,Y data for speed and for keeping zoom
% levels
% XX indicate on tuning curve which color is which
% XX add ability to switch from graded to distinguishable color scheme, also
% cyclical
% - add support for a second grouping variable - in that case, should use a
% different graded color scheme for each, overlay the traces in the tuning
% curve view, and also provide a 2-d image of tuning curve
% XX add support for plot labels (traces in psth, x-axis in tuning curve)
% XX option to change filter type
% XX expand calc window to avoid falling off at the edges
% XX option for error bars
% XX fix bug where the raster from the last unit stays when there are no new
% spikes to plot

fprintf(1, 'Controls:\n')
fprintf(1, '- left/right arrow: select previous/next cluster\n')
fprintf(1, '- up/down arrow: change smoothing of psth curves\n')
fprintf(1, '- c: dialog box to pick a new cluster ID number\n')
fprintf(1, '- t: toggle showing psth traces for each grouping variable or just the\n')
fprintf(1, 'overall. If showing just overall, raster is sorted chronologically. If\n')
fprintf(1, 'showing by grouping variable, raster is sorted by that variable.\n')
fprintf(1, '- r: select a new range within which to count spikes for the tuning curve\n')
fprintf(1, '- f: change smoothing filter type\n')
fprintf(1, '- x: change color scheme\n')
fprintf(1, '- e: turn on error bars\n')
fprintf(1, '- 1/2: decrease/increase the raster tick size\n')

if ~isempty(varargin)
    p = varargin{1}; % parameters supplied by user
else
    p = [];
end

if(plotmap)
    load('eyeCalibrated.mat')
    load('eyeRTs_continuous.mat')

myData.calibrated_eye_x = calibrated_eye_x;
myData.calibrated_eye_y = calibrated_eye_y;
myData.eyeRTs_continuous = eyeRTs_continuous;

end

myData.time_window = [-0.05, 0.05];

params.groupingName = getOr(p, 'groupingName', 'Grouping variable value');
params.groupingLegend = getOr(p, 'groupingLegend', unique(trGroups));
params.Z_Name = getOr(p, 'Z_Name', 'HI');
params.smoothSize = 15; % in msec, stdev of gaussian smoothing filter
params.clusterIndex = 1;
params.rasterScale = floor(numel(eventTimes)/100); % height of ticks
params.window = window;
params.showAllTraces = true;
params.showErrorShading = false;
params.startRange = window(1);
params.stopRange = window(2);
params.binSize = 0.001;
params.filterType = 1;
params.smWin = genSmWin(params);

spikeTimes= [];
clusters = [];
for i = 1:size(Z_Spikes.data, 2)

    if(~isempty(Z_Spikes.data(i).neural_timeStamps))
        spikeTimes = [spikeTimes, Z_Spikes.data(i).neural_timeStamps(1, :)];
        clusters = [clusters, Z_Spikes.data(i).neural_timeStamps(2, :)];
    end
end

[B, I] = sort(spikeTimes);
spikeTimes = spikeTimes(I);
myData.spikeTimes = spikeTimes;
clu = clusters(I);
myData.clu = clu;

if ~issorted(eventTimes)
    [eventTimes, ii] = sort(eventTimes);
    trGroups = trGroups(ii);
end

myData.eventTimes = eventTimes(:);
myData.trGroups = trGroups(:);
myData.clusterIDs = unique(clu);
myData.trGroupLabels = unique(myData.trGroups);
myData.nGroups = length(myData.trGroupLabels);
myData.plotAxes = [];
myData.plotclassification = plotclassification;
myData.plotdelayed = plotdelayed;
myData.plotmap = plotmap;
params.colorType = 1;
params.colors = genColors(params.colorType, 2);

myData.params = params;
myData.sp = Z_Spikes.sp;
myData.Z_Spikes = Z_Spikes;
myData.Z_Analog = Z_Analog;
f = figure; f.Color = 'w';
% set(f, 'Renderer', 'painters')

myData.flag = 0;
set(f, 'UserData', myData);
set(f, 'KeyPressFcn', @(f,k)psthViewerCallback(f, k));

% Create a text box
hTextBox = uicontrol('Style', 'edit', 'Position', [20 20 80 40], 'String', 'ROC');

% Create a button to save the text box input
hButton = uicontrol('Style', 'pushbutton', 'Position', [120 20 80 40], 'String', 'Save ROC', ...
    'Callback', @(src, event) saveTextCallback(src, event, f));

hTextBox2 = uicontrol('Style', 'edit', 'Position', [20 60 80 40], 'String', 'Unit Type');

hTextBox3 = uicontrol('Style', 'edit', 'Position', [20 100 80 40], 'String', 'RF Trusted');

hTextBox4 = uicontrol('Style', 'edit', 'Position', [20 140 80 40], 'String', 'win_st');

hTextBox5 = uicontrol('Style', 'edit', 'Position', [20 180 80 40], 'String', 'win_en');

% Create a button to save the text box input
hButton2 = uicontrol('Style', 'pushbutton', 'Position', [120 60 80 40], 'String', 'Save Unit Type', ...
    'Callback', @(src, event) saveText2Callback(src, event, f));

% Create a button to save the text box input
h1Button = uicontrol('Style', 'pushbutton', 'Position', [230 20 80 20], 'String', 'Save data', ...
    'Callback', @(src, event) exitCallback(src, event, f));

% Create a button to save the text box input
hButton4 = uicontrol('Style', 'pushbutton', 'Position', [120 140 80 40], 'String', 'set end', ...
    'Callback', @(src, event) saveText4Callback(src, event, f));

% Create a button to save the text box input
hButton5 = uicontrol('Style', 'pushbutton', 'Position', [120 180 80 40], 'String', 'set start', ...
    'Callback', @(src, event) saveText5Callback(src, event, f));


% Add text box data to the figure's UserData
myData.hTextBox = hTextBox;
myData.textBoxData = '';
set(f, 'UserData', myData);

% Add text box data to the figure's UserData
myData.hTextBox2 = hTextBox2;
myData.textBoxData2 = '';
set(f, 'UserData', myData);

% Add text box data to the figure's UserData
myData.hTextBox3 = hTextBox3;
myData.textBoxData3 = '';
set(f, 'UserData', myData);

% Add text box data to the figure's UserData
myData.hTextBox4 = hTextBox4;
myData.textBoxData4 = '';
set(f, 'UserData', myData);

% Add text box data to the figure's UserData
myData.hTextBox5 = hTextBox5;
myData.textBoxData5 = '';
set(f, 'UserData', myData);


% Make plots
if isempty(myData.plotAxes)
    for pidx = 1:10
        subplot(2,5,pidx);
        myData.plotAxes(pidx) = gca;
    end
    set(f, 'UserData', myData);
end
su = Z_Spikes.su;
if(~isfield(su, 'ROC'))
    for i = 1:size(su, 2)
        su(i).ROC = nan;
    end
end
if(~isfield(su, 'unitType'))
    for i = 1:size(su, 2)
        su(i).unitType = 0;
    end
end
myData.su = su;
su = myData.su;

set(f, 'UserData', myData);

myData.hTextBox3 = hTextBox3;
myData.textBoxData3 = '';
set(f, 'UserData', myData);

myData.sp = Z_Spikes.sp;
set(f, 'UserData', myData);


psthViewerPlot(f)

end

function psthViewerPlot(f)

% fprintf(1,'plot with fig %d\n', get(f,'Number'));
myData = get(f,'UserData');
p = myData.params;
sp = myData.sp;
if(isfield(myData.su, 'RF_Trusted') && myData.su(p.clusterIndex).RF_Trusted)
    set(myData.hTextBox3, 'String', 'Trusted');
else
    set(myData.hTextBox3, 'String', 'Not trusted');
end

if(~isfield(myData.su, 'Ipsi_contra') & myData.plotdelayed)
    myData.su(p.clusterIndex).Ipsi_contra = 0;
    set(myData.hTextBox3, 'String', 'contra');
elseif(isfield(myData.su, 'Ipsi_contra') & myData.plotdelayed & myData.su(p.clusterIndex).Ipsi_contra)
    set(myData.hTextBox3, 'String', 'ipsi');
    myData.su(p.clusterIndex).Ipsi_contra = 1;
else
    set(myData.hTextBox3, 'String', 'contra');
    myData.su(p.clusterIndex).Ipsi_contra = 0;
end
set(f, 'UserData', myData);

set(myData.hTextBox5, 'String', num2str(myData.time_window(1)*1000));
set(myData.hTextBox4, 'String', num2str(myData.time_window(2)*1000));

% pick the right spikes
st = myData.spikeTimes(myData.clu==myData.clusterIDs(myData.params.clusterIndex));

plotWindow = p.window;
calcWindow = plotWindow+p.smoothSize*3/1000*[0 1]; %half gaussian

tic
% compute everything
%[psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(st, myData.eventTimes, myData.params.window, 0.001);
[psth, bins, rasterX, rasterY, spikeCounts, ba] = psthAndBA(st, myData.eventTimes, calcWindow, p.binSize);

toc
trGroupLabels = myData.trGroupLabels;
nGroups = myData.nGroups;
inclRange = bins>p.startRange & bins<=p.stopRange;
spikeCounts = sum(ba(:,inclRange),2)./(p.stopRange-p.startRange);

% PSTH smoothing filter
smWin = p.smWin;

% smooth ba
baSm = conv2(smWin,1,ba', 'same')'./p.binSize;

% construct psth(s) and smooth it
if p.showAllTraces
    psthSm = zeros(nGroups, numel(bins));
    if p.showErrorShading
        stderr = zeros(nGroups, numel(bins));
    end
    for g = 1:nGroups
        theseTr = myData.trGroups==trGroupLabels(g);
        psthSm(g,:) = mean(baSm(theseTr,:));
        if p.showErrorShading
            stderr(g,:) = std(baSm(theseTr,:))./sqrt(sum(theseTr));
        end
    end
else
    psthSm = mean(baSm);
    if p.showErrorShading
        stderr = std(baSm)./sqrt(size(baSm,1));
    end

end

% compute raster
if p.showAllTraces
    [sortedGroups, inds] = sort(myData.trGroups);
    [tr,b] = find(ba(inds,:));
else
    [tr,b] = find(ba);
end
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything

% scale the raster ticks
rasterY(2:3:end) = rasterY(2:3:end)+(myData.params.rasterScale-1);

if isempty(rasterX); rasterX = NaN; rasterY = NaN; end % so that there's something to plot and the plot will clear

% compute the tuning curve
tuningCurve = zeros(nGroups,2);
for g = 1:nGroups
    theseCounts = spikeCounts(myData.trGroups==trGroupLabels(g));
    tuningCurve(g,1) = mean(theseCounts);
    tuningCurve(g,2) = std(theseCounts)./sqrt(length(theseCounts));
end

% ****************   Waveform over channels
myData = get(f, 'UserData');
subplot(2, 5, [1 6]); hold off
plot_waveform_overChannels(sp, 0, myData.clusterIDs(p.clusterIndex));
myData.flag = 1;
set(f, 'UserData', myData);

% colors = myData.params.colors;

colors = [0, 1, 0; 1, 0, 0; 0, 1, 0; 1, 0, 0; 0, 1, 0; 1, 0, 0; 0, 0, 1; 1, 0, 1; 0, 0, 1; 1, 0, 1; 0, 0, 1; 1, 0, 1];

for kk = 1:3
    % subplot(3,1,1);
    axes(myData.plotAxes(kk+1));
    hold off;
    if p.showAllTraces
        for g = (2*kk-1):(2*kk)
            if p.showErrorShading
                plotWithErr(bins, psthSm(g,:), stderr(g,:), colors(g,:));
            else
                plot(bins, psthSm(g,:), 'Color', colors(g,:), 'LineWidth', 2.0);
            end
            hold on;
        end
    else
        if p.showErrorShading
            plotWithErr(bins, psthSm, stderr, 'k');
        else
            plot(bins, psthSm, 'k');
        end
    end
    xlim(plotWindow);
    title(['Cluster ' num2str(myData.clusterIDs(p.clusterIndex)), '   Red: to the left (Ipsi)']);
    xlabel('Time (s)');
    ylabel('Firing rate (sp/s)');
    yl = ylim();
    hold on;
    plot(p.startRange*[1 1], yl, 'k--');
    plot(p.stopRange*[1 1], yl, 'k--');
    makepretty;
    box off;
    % subplot(3,1,2);
    axes(myData.plotAxes(kk+1));
    hold on;
    rasterscale = 150;
    if p.showAllTraces
        ug = unique(myData.trGroups);

        for g = (2*kk-1):(2*kk)
            theseTr = find(sortedGroups==ug(g));
            ry = rasterY>=theseTr(1) & rasterY<=theseTr(end);

            if any(ry)
                if(g == (2*kk-1))
                    ming = min(min(rasterY(ry)))/(rasterscale/max(max(psthSm((2*kk-1):(2*kk), :))));
                end
                ry(2:3:end) = ry(1:3:end);
                plot(rasterX(ry), rasterY(ry)/(rasterscale/max(max(psthSm((2*kk-1):(2*kk), :)))) - ming + max(max(psthSm((2*kk-1):(2*kk), :)))*1.1 , 'Color', colors(g,:));
            else
                plot(0, 0, 'Color', colors(g,:));
            end
            hold on;
        end
    else
        plot(rasterX,rasterY, 'k');
    end

    xlim(myData.params.window);
    xline(0, '--')
    ylim([0 max(rasterY(ry)/(rasterscale/max(max(psthSm((2*kk-1):(2*kk), :)))) - ming + max(max(psthSm((2*kk-1):(2*kk), :)))*1.1)]);
    ylabel('Event Number');
    xlabel('Time (s)');
    makepretty;
    box off;

end


% colors = [0, 0, 1; 1, 0, 1; 0, 0, 1; 1, 0, 1];
myData.plotclassification
if(myData.plotclassification)
    for kk = 4:5
        % subplot(3,1,1);
        axes(myData.plotAxes(kk+3));
        hold off;
        if p.showAllTraces
            for g = (2*kk-1):(2*kk)
                if p.showErrorShading
                    plotWithErr(bins, psthSm(g,:), stderr(g,:), colors(g,:));
                else
                    plot(bins, psthSm(g,:), 'Color', colors(g,:), 'LineWidth', 2.0);
                end
                hold on;
            end
        else
            if p.showErrorShading
                plotWithErr(bins, psthSm, stderr, 'k');
            else
                plot(bins, psthSm, 'k');
            end
        end
        xlim(plotWindow);
        title(['Cluster ' num2str(myData.clusterIDs(p.clusterIndex))]);
        xlabel('Time (s)');
        ylabel('Firing rate (sp/s)');
        yl = ylim();
        hold on;
        plot(p.startRange*[1 1], yl, 'k--');
        plot(p.stopRange*[1 1], yl, 'k--');
        makepretty;
        box off;
        % subplot(3,1,2);
        axes(myData.plotAxes(kk+3));
        hold on;
        rasterscale = 150;
        if p.showAllTraces
            ug = unique(myData.trGroups);

            for g = (2*kk-1):(2*kk)
                theseTr = find(sortedGroups==ug(g));
                ry = rasterY>=theseTr(1) & rasterY<=theseTr(end);

                if any(ry)
                    if(g == (2*kk-1))
                        ming = min(min(rasterY(ry)))/(rasterscale/max(max(psthSm((2*kk-1):(2*kk), :))));
                    end
                    ry(2:3:end) = ry(1:3:end);
                    plot(rasterX(ry), rasterY(ry)/(rasterscale/max(max(psthSm((2*kk-1):(2*kk), :)))) - ming + max(max(psthSm((2*kk-1):(2*kk), :)))*1.1 , 'Color', colors(g,:));
                else
                    plot(0, 0, 'Color', colors(g,:));
                end
                hold on;
            end
        else
            plot(rasterX,rasterY, 'k');
        end

        xlim(myData.params.window);
        xline(0, '--')
        ylim([0 max(rasterY(ry)/(rasterscale/max(max(psthSm((2*kk-1):(2*kk), :)))) - ming + max(max(psthSm((2*kk-1):(2*kk), :)))*1.1)]);
        ylabel('Event Number');
        xlabel('Time (s)');
        makepretty;
        box off;

        if(kk == 5)
            x = [myData.time_window(1), myData.time_window(1), myData.time_window(2), myData.time_window(2)];
            y = [0, max(rasterY(ry)/(rasterscale/max(max(psthSm((2*kk-1):(2*kk), :)))) - ming + max(max(psthSm((2*kk-1):(2*kk), :)))*1.1) max(rasterY(ry)/(rasterscale/max(max(psthSm((2*kk-1):(2*kk), :)))) - ming + max(max(psthSm((2*kk-1):(2*kk), :)))*1.1), 0];
            patch(x, y, 'red','FaceAlpha',.1)
        end

    end
elseif(myData.plotdelayed)
        for kk = 4:6
        % subplot(3,1,1);
        axes(myData.plotAxes(kk+3));
        hold off;
        if p.showAllTraces
            for g = (2*kk-1):(2*kk)
                if p.showErrorShading
                    plotWithErr(bins, psthSm(g,:), stderr(g,:), colors(g,:));
                else
                    plot(bins, psthSm(g,:), 'Color', colors(g,:), 'LineWidth', 2.0);
                end
                hold on;
            end
        else
            if p.showErrorShading
                plotWithErr(bins, psthSm, stderr, 'k');
            else
                plot(bins, psthSm, 'k');
            end
        end
        xlim(plotWindow);
        title(['Cluster ' num2str(myData.clusterIDs(p.clusterIndex))]);
        xlabel('Time (s)');
        ylabel('Firing rate (sp/s)');
        yl = ylim();
        hold on;
        plot(p.startRange*[1 1], yl, 'k--');
        plot(p.stopRange*[1 1], yl, 'k--');
        makepretty;
        box off;
        % subplot(3,1,2);
        axes(myData.plotAxes(kk+3));
        hold on;
        rasterscale = 150;
        if p.showAllTraces
            ug = unique(myData.trGroups);

            for g = (2*kk-1):(2*kk)
                theseTr = find(sortedGroups==ug(g));
                ry = rasterY>=theseTr(1) & rasterY<=theseTr(end);

                if any(ry)
                    if(g == (2*kk-1))
                        ming = min(min(rasterY(ry)))/(rasterscale/max(max(psthSm((2*kk-1):(2*kk), :))));
                    end
                    ry(2:3:end) = ry(1:3:end);
                    plot(rasterX(ry), rasterY(ry)/(rasterscale/max(max(psthSm((2*kk-1):(2*kk), :)))) - ming + max(max(psthSm((2*kk-1):(2*kk), :)))*1.1 , 'Color', colors(g,:));
                else
                    plot(0, 0, 'Color', colors(g,:));
                end
                hold on;
            end
        else
            plot(rasterX,rasterY, 'k');
        end

        xlim(myData.params.window);
        xline(0, '--')
        ylim([0 max(rasterY(ry)/(rasterscale/max(max(psthSm((2*kk-1):(2*kk), :)))) - ming + max(max(psthSm((2*kk-1):(2*kk), :)))*1.1)]);
        ylabel('Event Number');
        xlabel('Time (s)');
        makepretty;
        box off;

    end
end

% ROC
min_ROC_time = -100; max_ROC_time = 500;
ROC_time = min_ROC_time:max_ROC_time;
ROC = [];
time_re_target = (myData.params.window(1)*1000):(myData.params.window(2)*1000);

for i = 1:length(ROC_time)
    colmatch = find(time_re_target == ROC_time(i));
    % col1: time col2: AUC Value
    ROC(end+1,:) = [ROC_time(i), calcROC_BC(baSm(myData.trGroups==trGroupLabels(2),colmatch), baSm(myData.trGroups==trGroupLabels(1),colmatch))];
end
if(~myData.plotdelayed)
    axes(myData.plotAxes(9)); hold off;
else
    axes(myData.plotAxes(10)); hold off;
end
plot(ROC(:,1),ROC(:,2),'b-'); hold on;
title('ROC (red: detrended)');
axis([min_ROC_time max_ROC_time 0 1])
plot([0 0],[0 1],'k--')
plot([min_ROC_time max_ROC_time],[.5 .5],'k--')
plot([min_ROC_time max_ROC_time],[.6 .6],'k-.')
plot([min_ROC_time max_ROC_time],[.4 .4],'k-.')
onset_latency = find_onset_latency_ROC(ROC_time,ROC(:,2),[30 100],[.6 8 10],0);

if(~isnan(onset_latency))
    xline(onset_latency, '-.k')
    text(10,0.9,['ROC =',num2str(onset_latency)]);
end

baseline_ROC_start_index = find(ROC_time == -100);
baseline_ROC_end_index = find(ROC_time == 30);
if sum(ROC(baseline_ROC_start_index:baseline_ROC_end_index,2)) > 0
    % Then detrend the ROC from -100 to +30, replot, and
    % recalculate threshold
    R = linfit_jo(-100:30,ROC(baseline_ROC_start_index:baseline_ROC_end_index,2)');
    lineartrend = [];
    for x = min_ROC_time:max_ROC_time
        lineartrend(end+1) = R(1)*x+R(3)-0.5;   % to center it around 0.5
    end
    ROC_detrend = ROC(:,2)-lineartrend';

    if(~myData.plotdelayed)
        axes(myData.plotAxes(9)); hold on;
    else
        axes(myData.plotAxes(10)); hold on;
    end
    plot(ROC(:,1),ROC_detrend,'r-');
    %     title('Detrended time series ROC on spike density functions');
    axis([min_ROC_time max_ROC_time 0 1])
    plot([0 0],[0 1],'k--')
    plot([min_ROC_time max_ROC_time],[.5 .5],'k--')
    plot([min_ROC_time max_ROC_time],[.6 .6],'k-.')
    plot([min_ROC_time max_ROC_time],[.4 .4],'k-.')

    onset_latency_detrended = find_onset_latency_ROC(ROC_time,ROC_detrend,[30 100],[.6 8 10],0);
    if(~isnan(onset_latency_detrended))
        xline(onset_latency_detrended, '-.k')
        text(200,0.9,['det ROC =',num2str(onset_latency_detrended)]);
    end

    if ~isnan(onset_latency)
        onset_latency = min(onset_latency_detrended, onset_latency);
    end

end

curROC = myData.su(p.clusterIndex).ROC;
if(isnan(curROC))
    myData.su(p.clusterIndex).ROC = onset_latency;
end

set(myData.hTextBox2, 'String', num2str(myData.su(p.clusterIndex).unitType));
set(myData.hTextBox, 'String', num2str(myData.su(p.clusterIndex).ROC));
set(f, 'UserData', myData);

% ROC
min_ROC_time = -100; max_ROC_time = 500;
ROC_time = min_ROC_time:max_ROC_time;
ROC = [];
time_re_target = (myData.params.window(1)*1000):(myData.params.window(2)*1000);

for i = 1:length(ROC_time)
    colmatch = find(time_re_target == ROC_time(i));
    % col1: time col2: AUC Value
    ROC(end+1,:) = [ROC_time(i), calcROC_BC(baSm(myData.trGroups==trGroupLabels(1),colmatch), baSm(myData.trGroups==trGroupLabels(2),colmatch))];
end
onset_latency = find_onset_latency_ROC(ROC_time,ROC(:,2),[30 100],[.6 8 10],0);

if(~isnan(onset_latency))
    xline(onset_latency, '-.k')
    text(10,0.1,['ROC =',num2str(onset_latency)]);
end

baseline_ROC_start_index = find(ROC_time == -100);
baseline_ROC_end_index = find(ROC_time == 30);
if sum(ROC(baseline_ROC_start_index:baseline_ROC_end_index,2)) > 0
    % Then detrend the ROC from -100 to +30, replot, and
    % recalculate threshold
    R = linfit_jo(-100:30,ROC(baseline_ROC_start_index:baseline_ROC_end_index,2)');
    lineartrend = [];
    for x = min_ROC_time:max_ROC_time
        lineartrend(end+1) = R(1)*x+R(3)-0.5;   % to center it around 0.5
    end
    ROC_detrend = ROC(:,2)-lineartrend';

    %     axes(myData.plotAxes(6));hold on
    %     plot(ROC(:,1),ROC_detrend,'r-');
    % %     title('Detrended time series ROC on spike density functions');
    %     axis([min_ROC_time max_ROC_time 0 1])
    %     plot([0 0],[0 1],'k--')
    %     plot([min_ROC_time max_ROC_time],[.5 .5],'k--')
    %     plot([min_ROC_time max_ROC_time],[.6 .6],'k-.')
    %     plot([min_ROC_time max_ROC_time],[.4 .4],'k-.')

    onset_latency_detrended = find_onset_latency_ROC(ROC_time,ROC_detrend,[30 100],[.6 8 10],0);
    if(~isnan(onset_latency_detrended))
        xline(onset_latency_detrended, '-.k')
        text(200,0.1,['det ROC =',num2str(onset_latency_detrended)]);
    end

    if ~isnan(onset_latency)
        onset_latency = min(onset_latency_detrended, onset_latency);
    end

end

% myData = get(f, 'UserData');
% su = myData.su;

if(myData.plotmap)
    Z_Analog = myData.Z_Analog;
    Z_Spikes = myData.Z_Spikes;

    timestamp = (0:size(Z_Analog.data, 2)-1) ./ Z_Analog.info.SampleRate;

    spikeTimes= [];
    clusters = [];
    for i = 1:size(Z_Spikes.data, 2)
        if(~isempty(Z_Spikes.data(i).neural_timeStamps))
            spikeTimes = [spikeTimes, Z_Spikes.data(i).neural_timeStamps(1, :)];
            clusters = [clusters, Z_Spikes.data(i).neural_timeStamps(2, :)];
        end
    end

    cluster_want = myData.su(p.clusterIndex).clusterId + 1;
    calibrated_eye_x = myData.calibrated_eye_x;
    calibrated_eye_y = myData.calibrated_eye_y;

    % cluster_want = 35; % You can look at specific cluster map
    spikeTimess = spikeTimes(find(clusters == cluster_want));

    eyeRTs_continuous = myData.eyeRTs_continuous(2:end-2);
    spiking_timepoints = spikeTimess;  % spiking timepoints

    saccade_timepoints = timestamp(eyeRTs_continuous);  % saccade timepoints
    saccade_vectors = [calibrated_eye_x(eyeRTs_continuous + 70) - calibrated_eye_x(eyeRTs_continuous - 20); calibrated_eye_y(eyeRTs_continuous + 70) - calibrated_eye_y(eyeRTs_continuous - 20)]';

    % Remove saccade vectors where either x or y component is greater than 50
    sac_abs = sqrt(saccade_vectors(:, 1).^2 + saccade_vectors(:, 2) .^ 2);
    valid_saccades = ~(sac_abs > 70 | sac_abs < 0);
    saccade_timepoints = saccade_timepoints(valid_saccades);
    saccade_vectors = saccade_vectors(valid_saccades, :);

    % Define resolution
    resolution = 100;  % Number of bins for the heatmap

    % Bin the saccade vectors
    x_bins = linspace(-50, 50, resolution + 1);
    y_bins = linspace(-50, 50, resolution + 1);

    % Compute neural activity for each bin
    heatmap = zeros(resolution+1, resolution+1);
    counts = zeros(resolution+1, resolution+1);

    % Define the time window around each saccade onset to consider (e.g., 100 ms before to 100 ms after)
    time_window = myData.time_window;

    for i = 1:length(saccade_timepoints) % This is now aligned to target onset
        x = saccade_vectors(i, 1);
        y = saccade_vectors(i, 2);

        % Find the bin indices for the saccade vector
        [~, x_bin] = min(abs(x_bins - x));
        [~, y_bin] = min(abs(y_bins - y));

        % Find spikes in the time window around the saccade onset (first saccade after photodiode onset)
        onset = saccade_timepoints(i);
        %         onset = Z_Analog.PD_sec(i , 1);
        spikes_in_window = spiking_timepoints(spiking_timepoints >= (onset + time_window(1)) & spiking_timepoints <= (onset + time_window(2)));

        % Compute the neural activity (e.g., spike count)
        activity = length(spikes_in_window);

        % Update the heatmap and counts
        heatmap(y_bin, x_bin) = heatmap(y_bin, x_bin) + activity;
        counts(y_bin, x_bin) = counts(y_bin, x_bin) + 1;

    end

    % Normalize the heatmap by the counts
    heatmap = heatmap ./ counts;
    heatmap(isnan(heatmap)) = 0;  % Replace NaNs with 0

    % Apply Gaussian smoothing to the heatmap
    sigma = 1; % Standard deviation for Gaussian kernel
    smoothed_heatmap = imgaussfilt(heatmap, sigma);

    % Find the maximum bin
    [max_value, max_index] = max(smoothed_heatmap(:));
    [max_y, max_x] = ind2sub(size(smoothed_heatmap), max_index);
    max_x_coord = (x_bins(max_x) + x_bins(max_x+1)) / 2;
    max_y_coord = (y_bins(max_y) + y_bins(max_y+1)) / 2;

    % fprintf('Maximum bin coordinates: X = %.2f, Y = %.2f\n', max_x_coord*1.5, max_y_coord*1.5);

    % Method 2: Weighted Average of Top Bins
    num_top_bins = 5;  % Number of top bins to consider
    [sorted_values, sorted_indices] = sort(smoothed_heatmap(:), 'descend');
    top_indices = sorted_indices(1:num_top_bins);
    [top_y, top_x] = ind2sub(size(smoothed_heatmap), top_indices);

    top_x_coords = (x_bins(top_x) + x_bins(top_x+1)) / 2;
    top_y_coords = (y_bins(top_y) + y_bins(top_y+1)) / 2;
    weights = sorted_values(1:num_top_bins) / sum(sorted_values(1:num_top_bins));

    weighted_x_coord = sum(top_x_coords .* weights');
    weighted_y_coord = sum(top_y_coords .* weights');

    % Plot the heatmap
    subplot(2, 5, [5, 10]);hold off;
    imagesc(x_bins*1.5, y_bins*1.5, smoothed_heatmap);hold on
    set(gca, 'YDir', 'normal');
    % colormap('hot');
    colorbar;
    xlabel('Saccade Vector X');
    ylabel('Saccade Vector Y');
    title(['Coordinates: ', num2str(weighted_x_coord*1.5), ', ', num2str(weighted_y_coord*1.5)]);

    xline(weighted_x_coord*1.5)
    yline(weighted_y_coord*1.5)

    myData.Z_Spikes.su(p.clusterIndex).RF_x = weighted_x_coord*1.5;
    myData.Z_Spikes.su(p.clusterIndex).RF_y = weighted_y_coord*1.5;
end

set(f, 'UserData', myData);

end

function psthViewerCallback(f, keydata)

% fprintf('callback on %d with source %d\n', f.Number, keydata.Source.Number);

myData = get(f, 'UserData');
sp = myData.sp;
p = myData.params;

switch keydata.Key
    case 'rightarrow' % increment cluster index

        p.clusterIndex = p.clusterIndex+1;
        if p.clusterIndex>length(myData.clusterIDs)
            p.clusterIndex=1;
        end

    case 'leftarrow' % decrement cluster index

        p.clusterIndex = p.clusterIndex-1;
        if p.clusterIndex<1
            p.clusterIndex=length(myData.clusterIDs);
        end

    case 'uparrow' % increase smoothing
        p.smoothSize = p.smoothSize*1.2;
        p.smWin = genSmWin(p);

    case 'downarrow' % decrease smoothing
        p.smoothSize = p.smoothSize/1.2;
        p.smWin = genSmWin(p);

    case 'f'
        p.filterType = p.filterType+1;
        if p.filterType>3; p.filterType = 1; end

        p.smWin = genSmWin(p);

    case 't'
        if(myData.su(p.clusterIndex).RF_Trusted == 1)
            myData.su(p.clusterIndex).RF_Trusted = 0;
        else
            myData.su(p.clusterIndex).RF_Trusted = 1;
        end

    case 'i'
        if(myData.su(p.clusterIndex).Ipsi_contra == 1)
            myData.su(p.clusterIndex).Ipsi_contra = 0;
        else
            myData.su(p.clusterIndex).Ipsi_contra = 1;
        end

    case 'p'
        [myData.su(p.clusterIndex).RF_x, myData.su(p.clusterIndex).RF_y] = ginput(1);


        %     case 'e' % whether to show standard error as shading
        %         p.showErrorShading = ~p.showErrorShading;

        %     case 't' % whether to plot the psth trace for each condition or just the overall one
        %         p.showAllTraces = ~p.showAllTraces;



        %     case 'r'
        %         ax = subplot(3,1,1); title('click start and stop of range')
        %         %         [startRange,~] = ginput();
        %         %         [stopRange,~] = ginput();
        %         waitforbuttonpress;
        %         q = get(ax, 'CurrentPoint');
        %         p.startRange = q(1,1);
        %         waitforbuttonpress;
        %         q = get(ax, 'CurrentPoint');
        %         p.stopRange = q(1,1);
        %         if p.stopRange<p.startRange
        %             tmp = p.startRange;
        %             p.startRange = p.stopRange;
        %             p.stopRange = tmp;
        %         end

    case 'c'
        newC = inputdlg('cluster ID?');
        ind = find(myData.clusterIDs==str2num(newC{1}),1);
        if ~isempty(ind)
            p.clusterIndex = ind;
        end

        %     case 'x'
        %         p.colorType = p.colorType+1;
        %         p.colors = genColors(p.colorType, myData.nGroups);

        %     case '2' % increase tick size
        %         p.rasterScale = p.rasterScale*1.5;

        %     case '1' % decrease tick size
        %         p.rasterScale = p.rasterScale/1.5;
end

myData.params = p;
set(f, 'UserData', myData);

% plot with new settings
psthViewerPlot(f)


end



function smWin = genSmWin(p)

switch p.filterType
    case 1 % half gaussian, causal
        gw = gausswin(round(p.smoothSize*6*3),3);
        gw(1:round(numel(gw)/2)) = 0;
        smWin = gw./sum(gw);
        fprintf(1, 'filter is causal half-gaussian with stdev %.2f ms\n', p.smoothSize*3);
    case 2 % gaussian
        gw = gausswin(round(p.smoothSize*6),3);
        smWin = gw./sum(gw);
        fprintf(1, 'filter is gaussian with stdev %.2f ms\n', p.smoothSize);
    case 3 % box
        smWin = ones(round(p.smoothSize*3),1);
        fprintf(1, 'filter is box with width %.2f ms\n', p.smoothSize*3);
end
end

function makepretty()
% set some graphical attributes of the current axis

set(get(gca, 'XLabel'), 'FontSize', 15);
set(get(gca, 'YLabel'), 'FontSize', 15);
set(gca, 'FontSize', 13);

set(get(gca, 'Title'), 'FontSize', 15);

ch = get(gca, 'Children');

for c = 1:length(ch)
    thisChild = ch(c);
    if strcmp('line', get(thisChild, 'Type'))
        if strcmp('.', get(thisChild, 'Marker'))
            set(thisChild, 'MarkerSize', 15);
        end
        if strcmp('-', get(thisChild, 'LineStyle'))
            set(thisChild, 'LineWidth', 2.0);
        end
    end
end

end

function colors = genColors(colorType, nColors)

switch mod(colorType, 3)
    case 1 % blue linear map
        colors = copper(nColors);
        colors = colors(:, [3 2 1]);
    case 2 % distinguishable colors - default matlab order with black and gray
        colors = get(gca, 'ColorOrder');
        colors = [colors; 0 0 0; 0.5 0.5 0.5];
        nc = size(colors,1);
        colors = colors(mod(0:nColors-1,nc)+1,:);
    case 0 % cyclical
        %colors = hsv(nColors);
        m = colorcet('C6');
        colors = zeros(nColors,3);
        for c = 1:3
            qidx = linspace(1,size(m,1), nColors+1);
            colors(:,c) = interp1(1:size(m,1), m(:,c), qidx(1:end-1));
        end
end
end

% Callback function to save the text box input
function saveTextCallback(src, event, f)
myData = get(f, 'UserData');
p = myData.params;
myData.textBoxData = get(myData.hTextBox, 'String');
myData.su(p.clusterIndex).ROC = str2num(get(myData.hTextBox, 'String'));
set(f, 'UserData', myData);
%     disp(['Saved text: ', myData.textBoxData]); % Display the saved text in the command window
end

% Callback function to save the text box input
function saveText2Callback(src, event, f)
myData = get(f, 'UserData');
p = myData.params;
myData.textBoxData2 = get(myData.hTextBox2, 'String');
myData.su(p.clusterIndex).unitType = str2num(get(myData.hTextBox2, 'String'));
set(f, 'UserData', myData);
end

% Callback function to save the text box input
function saveText4Callback(src, event, f)
myData = get(f, 'UserData');
p = myData.params;
myData.textBoxData4 = get(myData.hTextBox4, 'String');
myData.time_window(2) = str2num(get(myData.hTextBox4, 'String')) / 1000;
set(f, 'UserData', myData);
psthViewerPlot(f)

end

% Callback function to save the text box input
function saveText5Callback(src, event, f)
myData = get(f, 'UserData');
p = myData.params;
myData.textBoxData5 = get(myData.hTextBox5, 'String');
myData.time_window(1) = str2num(get(myData.hTextBox5, 'String')) / 1000;
set(f, 'UserData', myData);
psthViewerPlot(f)

end


% Callback function to save the text box input
function exitCallback(src, event, f)
myData = get(f, 'UserData');
su = myData.su;
Z_Name = myData.params.Z_Name;
myData.Z_Spikes.su = su;
Z_Spikes = myData.Z_Spikes;
save(Z_Name, 'Z_Spikes')
disp('Saved')
end