function do_posthoc(Folders, meta)

sesName = meta.sesName;

% set paths for where to find your data

myKsDir = Folders.datFolder;

eventFolder = [Folders.KiloFolder, sesName, '\'];

load([eventFolder, 'event_times_r.mat'])
load([eventFolder, 'event_times_l.mat'])

event_times_r = event_times_r{5}(1, :);
event_times_l = event_times_l{5}(1, :);

% Loading data from kilosort/phy easily

spp = loadKSdir(myKsDir);

% Plotting a driftmap

[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths, 'mark');
figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths, 'show');

% basic quantification of spiking plot

depthBins = 0:40:5000;
ampBins = 0:30:min(max(spikeAmps),800);
recordingDur = spp.st(end);

[pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);
sgtitle('Spike Amps over depth')

% Computing some useful details about spikes/neurons (like depths)

[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(spp.temps, spp.winv, spp.ycoords, spp.spikeTemplates, spp.tempScalingAmps);

% Looking at PSTHs aligned to some event

% if you now have a vector of relevant event times, called eventTimes (but
% not the cell array as above, just a vector):

window = [-0.3 1]; % look at spike times from 0.3 sec before each event to 1 sec after

% if your events come in different types, like different orientations of a
% visual stimulus, then you can provide those values as "trial groups",
% which will be used to construct a tuning curve. Here we just give a
% vector of all ones.
event_times = [event_times_r, event_times_l];

trialGroups = [zeros(size(event_times_r))+2, zeros(size(event_times_l))+3];

psthViewer(spp.st, spp.clu, event_times, window, trialGroups);

% use left/right arrows to page through the clusters


% PSTHs across depth

depthBinSize = 150; % in units of the channel coordinates, in this case µm
timeBinSize = 0.001; % seconds
bslWin = [-0.2 -0.05]; % window in which to compute "baseline" rates for normalization
psthType = 'norm'; % show the normalized version
eventName = 'stimulus onset'; % for figure labeling

[timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimes, spikeDepths, ...
    depthBinSize, timeBinSize, event_times_r, window, bslWin);

figure;
plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);

end