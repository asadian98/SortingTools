% Example script for some of the functions in the spikes repository. 
%
% These functions make it easy to work with spiking data from kilosort/phy 
% (or, for many functions, from anywhere)
%
% Please make an issue if this function is missing some dependencies.
%
% See https://github.com/cortex-lab/neuropixels/wiki/Other_analysis_methods
% for more explanation. 

%% add the repositories to your path

addpath(genpath('C:\...\github\spikes'))
addpath(genpath('C:\...\github\npy-matlab'))


%% set paths for where to find your data

myKsDir = 'C:\Users\CorneilLab\Desktop\AA\Kilosort\Belle\Be240402\dat';

myEventTimes = load('C:\Users\CorneilLab\Desktop\AA\Kilosort\Belle\Be240402\event_times_r.mat'); % a vector of times in seconds of some event to align to

load('C:\Users\CorneilLab\Desktop\AA\Kilosort\Belle\Be240402\event_times_r.mat')
load('C:\Users\CorneilLab\Desktop\AA\Kilosort\Belle\Be240402\event_times_l.mat')
%% Loading data from kilosort/phy easily

sp = loadKSdir(myKsDir)

% sp.st are spike times in seconds
% sp.clu are cluster identities
% spikes from clusters labeled "noise" have already been omitted

%% Plotting a driftmap

[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths, 'mark');
figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths, 'show');

%% basic quantification of spiking plot

depthBins = 0:40:4500;
ampBins = 0:30:min(max(spikeAmps),800);
recordingDur = sp.st(end);

[pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);

%% Computing some useful details about spikes/neurons (like depths)

[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);


%% Looking at PSTHs aligned to some event

% if you now have a vector of relevant event times, called eventTimes (but
% not the cell array as above, just a vector):

window = [-0.3 1]; % look at spike times from 0.3 sec before each event to 1 sec after

% if your events come in different types, like different orientations of a
% visual stimulus, then you can provide those values as "trial groups",
% which will be used to construct a tuning curve. Here we just give a
% vector of all ones. 
event_times = [event_times_r, event_times_l];

trialGroups = [zeros(size(event_times_r))+2, zeros(size(event_times_l))+3]; 

psthViewer(sp.st, sp.clu, event_times, window, trialGroups);

% use left/right arrows to page through the clusters


%% PSTHs across depth

depthBinSize = 150; % in units of the channel coordinates, in this case µm
timeBinSize = 0.001; % seconds
bslWin = [-0.2 -0.05]; % window in which to compute "baseline" rates for normalization
psthType = 'norm'; % show the normalized version
eventName = 'stimulus onset'; % for figure labeling

[timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimes, spikeDepths, ...
    depthBinSize, timeBinSize, event_times_r, window, bslWin);

figure;
plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);


%% Loading raw waveforms

% To get the true waveforms of the spikes (not just kilosort's template
% shapes), use the getWaveForms function:

gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
apD = dir(fullfile(myKsDir, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = 'Be240402_001.dat';         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes = ceil(sp.st*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu;

wf = getWaveForms(gwfparams);

figure; 
imagesc(squeeze(wf.waveFormsMean(2, :, :)))
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number'); 
colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;
