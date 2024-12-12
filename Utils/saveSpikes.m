function saveSpikes(Folders, meta, Physio_Info, saveWave_flag)

sesName = meta.sesName;
region = meta.region;
fileName = meta.fileName;
nCh = meta.nCh;

disp('Getting Sorted spikes ...')

getPhydata(sesName, Folders.sortingFolder, region)

cd(Folders.sortingFolder)
load([sesName, '_Phy_', region])

for i = 1:Physio_Info.numCh
    spk_arr = [];
    ID_arr = [];
    ch_cluster = find(cluster_channel == i);
    for k = 1:length(ch_cluster)
        spk_arr = [spk_arr, spikes_phy(2, find(spikes_phy(1, :) == (clusterId(1, ch_cluster(k)))))/meta.fs];
        ID_arr = [ID_arr, spikes_phy(1, find(spikes_phy(1, :) == (clusterId(1, ch_cluster(k)))))];
    end

    [spk_arr, I] = sort(spk_arr);
    Z_Spikes.data(i).neural_timeStamps = [spk_arr; ID_arr(I)];
    Z_Spikes.data(i).Units = unique(ID_arr(I));
    Z_Spikes.data(i).nUnits = length(unique(ID_arr(I)));
    Z_Spikes.data(i).scores = clusterScore(find(cluster_channel == i));
end
Z_Spikes.info.SampleRate = meta.fs;

% Extract waveforms and quality info and stuff ...
if(saveWave_flag)
    [sp, wf, wfOnPeak] = getSp(Folders.sortingFolder, Folders.save_dir, sesName, [fileName, '.dat'], nCh, 'waves', true);
else
    [sp, wf, wfOnPeak] = getSp(Folders.sortingFolder, Folders.save_dir, sesName, [fileName, '.dat'], nCh, 'waves', false);
end

sp2clust(sp, Folders.sortingFolder)

% make sure there's a folder:
if ~exist(fullfile(Folders.sortingFolder, 'figures'), 'dir')
    mkdir(fullfile(Folders.sortingFolder, 'figures'));
end

% set options:
opts.saveFigs = 1;
opts.dirFigs = Folders.sortingFolder;

% make dem figures and save:
%                 mkfig.waveform_overChannels(sp, opts);
%                 mkfig.waveformOverChannels_perCluster(sp, opts);
%                 mkfig.waveformAndSpikeCount_overChannels(sp, opts)
su = sp2su(sp, Folders.sortingFolder);

disp('Saving Z_Spikes ...')
cd(Folders.save_dir)
Z_Spikes.su = su;
Z_Spikes.sp = sp;

save(['Z_', sesName, '_', region, '_Spikes'], 'Z_Spikes','-v7.3')

Z_Wave.data = wf;

cd(Folders.save_dir)
if(saveWave_flag)
    disp('Saving Z_Wave ... ')
    save(['Z_', sesName, '_', region, '_Wave'], 'Z_Wave','-v7.3')
end
%         mkfig.unitSummary(su, opts);


% % Plot first 1000 waveforms for the specified unitID
% figure
% wf_time = linspace(-0.5, 1.25, length(-round(0.5/1000*30000):round(1.25/1000*30000)));
% UnitID = 55;
% UnitIdx = find([su.clusterId] == UnitID);
% for i = 1:100
% plot(wf_time, squeeze(wf{1, UnitIdx}(i, su(UnitIdx).peakCh, :) - mean(wf{1, UnitIdx}(i, su(UnitIdx).peakCh, :))), 'Color', [.7 .7 .7])
% hold on
% end
% plot(wf_time, su(UnitIdx).medWfOnPeakCh - mean(su(UnitIdx).medWfOnPeakCh), 'k')
% xlabel('time (ms)')
% ylabel('Amp (uV)')
% title(['Cluster ID: ', num2str(UnitID)])

end