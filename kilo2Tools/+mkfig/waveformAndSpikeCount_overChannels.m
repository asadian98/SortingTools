function hFig = waveformAndSpikeCount_overChannels(sp, opts)
%   hFig = mkfig.waveformAndSpikeCount_overChannels(sp, opts)

%%
if ~exist('opts', 'var')
    opts = [];
end

%%
hFig = figure('Position', [100 100 300 900]); 
figSz = [3 8];

%%
subplot(121); 
plot_waveform_overChannels(sp);

%%
subplot(122); 
plot_spikeCount_overChannles(sp);

%%

sgtitle(sp.info.dsn)
formatFig(hFig, figSz, 'default')
if isfield(opts, 'saveFigs') && opts.saveFigs == true
    if ~isfield(opts, 'dirFigs')
        opts.dirFigs = pwd;
    end
    saveas(hFig, fullfile(opts.dirFigs, 'figures', 'waveformAndSpikeCount_overChannels.pdf'));
end



