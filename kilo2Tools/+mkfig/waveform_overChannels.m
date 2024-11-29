function hFig = waveform_overChannels(sp, opts)
%   hFig = mkfig.waveform_overChannels(sp, opts)

if ~exist('opts', 'var')
    opts = [];
end

%%
hFig = figure('Position', [100 100 300 900]); 
figSz = [3 8];

plot_waveform_overChannels(sp);

%%
% sgtitle(sp.info.dsn)
formatFig(hFig, figSz)
if isfield(opts, 'saveFigs') && opts.saveFigs == true
    if ~isfield(opts, 'dirFigs')
        opts.dirFigs = pwd;
    end
    saveas(hFig, fullfile(opts.dirFigs, 'figures', 'waveform_overChannels.pdf'));
end
