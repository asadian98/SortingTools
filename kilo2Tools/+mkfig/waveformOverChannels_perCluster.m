function hFig = waveformOverChannels_perCluster(sp, opts)
%   hFig = mkfig.waveformOverChannels_perCluster(sp, opts)
%
% makes figure with a subplot per cluster. each subplot has the median
% waveform on each of the nChannels. works on kilosort output sp struct
%
% INPUT:
%   sp  - works on kilosort output sp struct (see getSp.m)
%%

if ~exist('opts', 'var')
    opts = [];
end

%%
hFig = figure('Position', [100 100 300 900]); 
figSz = [length(find(sp.clusterScore == 2)), 8];

yScaler = 1.2;
x  = 1:size(sp.medWfs, 3);

nX = numel(unique(sp.xcoords));
nY = numel(sp.ycoords);

newarray = find(sp.clusterScore == 2);

[~, idxSort] = sort(sp.peakCh(newarray)); % Just plot good clusters

iPlot = 1;
for iClu = newarray(idxSort)
    subplot(1, length(find(sp.clusterScore == 2)), iPlot)
    hold on
    for iCh = 1:nY
%         plot(((iCh-1)*yBuffer) + squeeze(sp.medWFs(iClu, iCh,:))');
        plot(x + sp.xcoords(iCh) - 0.5*length(x), sp.ycoords(iCh) + yScaler * squeeze(sp.medWfs(iClu, iCh,:))');
    end
    xlim([min(sp.xcoords) - 50, max(sp.xcoords) + 50])
    ylim([min(sp.ycoords) - 50, max(sp.ycoords) + 50])
    set(gca, 'XTick',[], 'YTick',[])
    if iPlot==1
        ylabel('y distance (µ)')
        xlabel('x distance (µ)')
        set(gca, 'XTick', unique(sp.xcoords), 'YTick', unique(sp.ycoords))
    end
    title(['cid: ' num2str(sp.clusterId(iClu))])
    iPlot = iPlot + 1;
end

%%
formatFig(hFig, figSz, 'default');
% nature, jnsci, eneuro, poster, default
sgtitle(sp.info.dsn)
if isfield(opts, 'saveFigs') && opts.saveFigs == true
    if ~isfield(opts, 'dirFigs')
        opts.dirFigs = pwd;
    end
    saveas(hFig, fullfile(opts.dirFigs, 'figures', 'waveformOverChannels_perCluster.pdf'));
end

