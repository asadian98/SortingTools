%%
% if you wish to get waveforms (time consuming), this is where it happens:

clc; clear;
datapath = 'C:\Users\CorneilLab\Desktop\AA\Belle';
ses_list = {'Be231124', 'Be231128', 'Be231129', 'Be231130', 'Be231201', ...
    'Be231205', 'Be231206', 'Be231207', 'Be231208', 'Be231211', 'Be231212', 'Be231213', 'Be240510', 'Be240514'};
regions = {'SC'};

folder_path = 'C:\Users\CorneilLab\Desktop\SortingTools';
addpath(genpath(folder_path));
folder_path = 'C:\Users\CorneilLab\Desktop\SortingTools\bombcell-main';
addpath(genpath(folder_path));

for ses_list_idx = 11:length(ses_list)
    ses_list_idx

    for j = 1:length(regions)

        cd(datapath)
        cd(ses_list{ses_list_idx})

        if(exist(['Z_', ses_list{ses_list_idx}, '_', regions{j}, '_Spikes.mat']))
            load(['Z_', ses_list{ses_list_idx}, '_', regions{j}, '_Spikes'])

            ksDir = [datapath, '\', ses_list{ses_list_idx}, '\dat_SC'];
            sp = Z_Spikes.sp;

            wf = [];
            clu = sp.spikeClusters;
            cids = unique(clu);
            nClu = length(cids);
            wfOnPeak = cell(1, nClu);
            p.Results.waves = 1;

            if p.Results.waves

                disp('Retreiving all (yes all) waveforms. This might take a while...')

                st  = sp.spikeTimesSamps;
                Fs = sp.sample_rate;
                datPath = fullfile(ksDir, sp.dat_path);
                dataType = sp.dtype;
                d = dir([ksDir '/*.dat']);
                nSamp = d.bytes/2/sp.n_channels_dat;
                dataSize = [sp.n_channels_dat nSamp];
                chanMap = readNPY(fullfile(ksDir, 'channel_map.npy'));
                %     gain = 0.6/512/500*1e6; % raw file units to uV ***SHOULD BE RIG SPECIFIC. NEED TO DO THIS...
                method = 'samples';

                % window is -0.5 to 1.25ms
                % -0.5:1
                wfWin = -round(0.5/300*Fs):round(1/300*Fs);
                nWFsamps = numel(wfWin);

                nChInFile = dataSize(1);
                nSamp = dataSize(2);
                mmf = memmapfile(datPath, 'Format', {dataType, [nChInFile nSamp], 'x'});

                nCh = length(chanMap);
                hWait = waitbar(0, 'Extracting individual waveforms....');
                wf = cell(1, nClu);
                tic
                for iClu = 1:nClu
                    fprintf(1, 'cluster %d (%d/%d)\n', cids(iClu), iClu, nClu);
                    theseST = st(clu==cids(iClu));
                    nWFsToLoad = length(theseST);
                    switch method
                        case 'secs'
                            extractST = round(theseST(randperm(length(theseST), nWFsToLoad))*Fs);
                        case 'samples'
                            extractST = round(theseST(randperm(length(theseST), nWFsToLoad)));
                    end
                    % in case a spike is detected at the very start or end of the
                    % recording its full waveform can't be extracted, so remove it:
                    extractST((extractST - wfWin(1))< 1) = [];
                    extractST((extractST + wfWin(end))> size(mmf.Data.x,2)) = [];
                    nWf = numel(extractST);


                    % get the waveforms:
                    theseWF = zeros(nWf, nCh, nWFsamps);
                    for i=1:nWf
                        tempWF = mmf.Data.x(1:nChInFile, extractST(i) + wfWin(1):extractST(i) + wfWin(end));
                        theseWF(i,:,:) = tempWF(chanMap+1,:);
                    end

                    wf{iClu}  = theseWF;
                    waitbar(iClu/nClu, hWait);
                    disp(toc)
                end
                %     save(fullfile(ksDir, 'sp_waveforms.mat'), 'wf', '-v7.3');

                % now extract only waveforms from the peak channel, and save:
                for iClu = 1:nClu
                    wfOnPeak{iClu} = squeeze(wf{iClu}(:, find(chanMap+1 == sp.peakCh(iClu)), :));
                end
                save(fullfile(ksDir, 'sp_waveformsOnPeak.mat'), 'wfOnPeak', '-v7.3');
                close(hWait)
            end
        end
    end
end



% % Code for computing ISI violations
% iClu = 8;
% for iClu = 15
% refDur = 0.001;     % estiamtion of refractory period duration (lnk: was 0.0015)
% minISI = 0.0005;    % min possible ISI given waveform length (lnk: this is only used to compute the false positive rate);
% [fpRate, numViolations] = ISIViolations(sp.spikeTimesSecs(sp.spikeClusters==cids(iClu)), minISI, refDur);
% end
%
% (numViolations/length(sp.spikeTimesSecs(sp.spikeClusters==cids(iClu))))*100
%
% isis = diff(sp.spikeTimesSecs(sp.spikeClusters==cids(iClu)));
% violations = isis <= 0.0015;
% num_violations = sum(violations);
% total_isis = length(isis);
% percentage_violations = (num_violations / total_isis) * 100;
%
% % Display results
% fprintf('Number of ISI violations: %d\n', num_violations);
% fprintf('Total number of ISIs: %d\n', total_isis);
% fprintf('Percentage of ISI violations: %.2f%%\n', percentage_violations);
%
% % Plot the ISI distribution
% figure;
% histogram(isis, 'BinWidth', 0.002, 'Normalization', 'probability'); % Adjust bin width as needed
% xlabel('Inter-Spike Interval (ms)');
% ylabel('Probability');
% title('ISI Distribution');
% grid on;

% Code for plotting waveforms on the peak channel
% Create a figure
% figure;
% hold on;
% data = wfOnPeak{iClu};
% % Define spacing between vertical axes
% vertical_spacing = 20;
% 
% % Loop through the second dimension (M)
% for m = 1
%     % Offset for vertical separation
%     vertical_offset = (m - 1) * vertical_spacing;
% 
%     % Loop through the first dimension (N)
%     for n = 1:1000
%         % Plot the waveform with vertical offset
%         plot(linspace(-0.5, 1, 151), movmean(wfOnPeak{iClu}(n, :) - mean(wfOnPeak{iClu}(n, :)), 3), 'color', [0.7, 0.7, 0.7]);
%         hold on
%     end
%     plot(linspace(-0.5, 1, 151), median(wfOnPeak{iClu} - mean(wfOnPeak{iClu}, 2)), 'color', [1, 0, 0]);
% end
% 
% % Add labels and legend
% xlabel('Samples (3rd dimension)');
% ylabel('Amplitude (with vertical offsets)');
% title('3D Matrix Visualization');
% hold off;


%%
% This script is for classification of cells based on their waveform and
% ISI

% First, find good units and put the average firing rate in a matrix



% I need to align the waveforms as well



clc; clear;
datapath = 'C:\Users\CorneilLab\Desktop\AA\Belle';
ses_list = {'Be231124', 'Be231128', 'Be231129', 'Be231130', 'Be231201', ...
    'Be231205', 'Be231206', 'Be231207', 'Be231208', 'Be231211', 'Be231212', 'Be231213', 'Be240510', 'Be240514'};
regions = {'SC'};
% datapath = 'C:\Users\CorneilLab\Desktop\AA\Grover';
% ses_list = {'06242022', '07072022', '07112022', '07122022', '07132022', ...
%     '07142022', '07152022', '07192022', '07212022', 'Gr230321_002', 'Gr230322', 'Gr230324', 'Gr230328', 'Gr230329', 'Gr230330', 'Gr230331', 'Gr230403', 'Gr230405'};
% regions = {'SC'};

% Belle cortex
% datapath = 'C:\Users\CorneilLab\Desktop\AA\Belle';
% ses_list = {'Be240229', 'Be240305', 'Be240306', 'Be240307', 'Be240308', ...
%     'Be240327', 'Be240328', 'Be240402', 'Be240403', 'Be240404', 'Be240405', 'Be240417', 'Be240418', 'Be240419', 'Be240423', 'Be240424', 'Be240426'};
% regions = {'PMd', 'M1'};

folder_path = 'C:\Users\CorneilLab\Desktop\SortingTools';
addpath(genpath(folder_path));
folder_path = 'C:\Users\CorneilLab\Desktop\SortingTools\bombcell-main';
addpath(genpath(folder_path));


waveformsMed_all = zeros(1, 151);
spikes = {};
depth = [];
unitType = [];
k = 1;
for ses_list_idx = 1:length(ses_list)
    ses_list_idx
    cd(datapath)
    cd(ses_list{ses_list_idx})

    for j = 1:length(regions)

        if(exist(['Z_', ses_list{ses_list_idx}, '_', regions{j}, '_Spikes.mat']))

            load(['Z_', ses_list{ses_list_idx}, '_', regions{j}, '_Spikes'])
            cd(['dat_', regions{j}])
            load('sp_waveformsOnPeak.mat')

            goodUnits = find([Z_Spikes.su(:).clusterScore] == 2);

            waveformsMed = zeros(1, 151); % 151 is fixed, number of timepoints for a given waveform
            for i = 1:length(goodUnits)
                waveformsMed(i, :) = median(wfOnPeak{goodUnits(i)} - mean(wfOnPeak{goodUnits(i)}, 2));
            end

            waveformsMed_all = [waveformsMed_all; waveformsMed];

            % Get all spikes for this session
            spikes_all = [];

            for ch = 1:size(Z_Spikes.data, 2)
                spikes_all = [spikes_all, Z_Spikes.data(ch).neural_timeStamps];
            end

            for i = 1:length(goodUnits)

                spikes{k} = spikes_all(1, find(spikes_all(2, :) == Z_Spikes.su(goodUnits(i)).clusterId+1));

                for ch = 1:size(Z_Spikes.data, 2)
                    if(ismember(Z_Spikes.su(goodUnits(i)).clusterId+1, Z_Spikes.data(ch).Units))
                        if(isfield(Z_Spikes.info, 'surface'))
                            depth(k) = (ch - Z_Spikes.info.surface) * 150;
                        else
                            depth(k) = ch*150;
                        end
                        break;
                    end
                end

                if(isfield(Z_Spikes.su(1), 'unitType'))
                    unitType(k) = Z_Spikes.su(goodUnits(i)).unitType;
                else
                    unitType(k) = 0;
                end

                k = k + 1;
            end

        end

    end
end

% Some waveforms are all zero --> remove them
waveformsMed_all(find(all(waveformsMed_all == 0, 2)), :) = [];

waveformsMed_all_norm = (waveformsMed_all - min(waveformsMed_all, [], 2)) ./ (max(waveformsMed_all, [], 2) - min(waveformsMed_all, [], 2));

% remove non-somatic waveforms
figure
non_somatic = zeros(1, size(waveformsMed_all, 1));
for non_somatic_idx = 1:size(waveformsMed_all, 1)
    waveform = waveformsMed_all(non_somatic_idx, :);

    % Original params from bombcell
    param.minThreshDetectPeaksTroughs = 0.2;
    param.minMainPeakToTroughRatio = 5;
    param.minWidthFirstPeak = 4;
    param.minWidthMainTrough = 5;
    param.firstPeakRatio = 3;
    param.minTroughToPeakRatio = 0.8;

    [nPeaks, nTroughs, mainPeak_before_size, mainPeak_after_size, mainTrough_size, ...
        mainPeak_before_width, width_after, mainTrough_width, peakLocs, troughLocs, PKS, TRS, troughLoc] = bc.qm.helpers.getWaveformPeakProperties(waveform, param);

    isNonSomatic = (abs(mainTrough_size ./ mainPeak_before_size) < param.minMainPeakToTroughRatio &...
        mainPeak_before_width < param.minWidthFirstPeak &...
        mainTrough_width < param.minWidthMainTrough &...
        abs(mainPeak_before_size./mainPeak_after_size) > param.firstPeakRatio) | ...
        abs(max([mainPeak_before_size, mainPeak_after_size], [], 2)./ mainTrough_size) > param.minTroughToPeakRatio;
    if(isNonSomatic)
        plot(waveform, 'r'); hold on;
    else
        plot(waveform, 'b'); hold on;
    end
    non_somatic(non_somatic_idx) = isNonSomatic;
end
title('Blue: Somatic,    Red: Non-somatic')


waveformsMed_all2 = waveformsMed_all(find(~non_somatic), :);

[N, m] = size(waveformsMed_all2); % Get the size of the input matrix
target_index = ceil(m / 3); % Define the target index to align the troughs to
aligned_waveforms = zeros(N, m); % Initialize the aligned waveform matrix

for i = 1:N
    waveform = movmean(waveformsMed_all2(i, :), 3); % Extract the waveform
    [~, trough_index] = min(waveform); % Find the index of the trough
    shift_amount = target_index - trough_index; % Calculate the shift amount

    waveform = movmean(waveformsMed_all2(i, :), 3);

    % Shift the waveform
    if shift_amount > 0
        % Shift right and pad with zeros at the beginning
        aligned_waveforms(i, :) = [zeros(1, shift_amount), waveform(1:end-shift_amount)];
    elseif shift_amount < 0
        % Shift left and pad with zeros at the end
        aligned_waveforms(i, :) = [waveform(-shift_amount+1:end), zeros(1, -shift_amount)];
    else
        % No shift needed
        aligned_waveforms(i, :) = waveform;
    end
end

waveformsMed_all_norm = (aligned_waveforms - min(aligned_waveforms, [], 2)) ./ (max(aligned_waveforms, [], 2) - min(aligned_waveforms, [], 2));

plot(waveformsMed_all_norm')

% % align spikes in time via depolarization trough
% [b, a] = butter(4, 250 / (1000 / 2), 'high');
% waveformsMed_all_filt = zeros(size(waveformsMed_all));
% for i = 1:size(waveformsMed_all, 1)
%     waveformsMed_all_filt(i, :) = filter(b, a, waveformsMed_all(i, :));
% end
% x = linspace(-0.5, 1.25, 54);
% waveformsMed_all_filt = movmean(waveformsMed_all, 1, 1);
% plot(x, waveformsMed_all(10, :))
% hold on
% plot(x, waveformsMed_all_filt(10, :))

log10isi_all = zeros(length(spikes), 101);
for i = 1:length(spikes)

    ISI = diff(sort(spikes{i}));
    log10isi = log10(ISI);
    x_range = -3:0.04:1;
    [f, xi] = ksdensity(log10isi, x_range);
    log10isi_all(i, :) = f;

end
% intervals = -3:0.04:1;
% intervals2 = intervals(1:end-1)+.02;
% log10_bins = 10.^intervals2';
% log10isi = log10isi.*(10.^intervals2)';
%%

waveforms_somatic = waveformsMed_all_norm;
waveforms_nonsomatic = waveformsMed_all(find(non_somatic), :);

spikes_somatic = spikes(find(~non_somatic));
spikes_nonsomatic = spikes(find(non_somatic));

log10isi_somatic = log10isi_all(find(~non_somatic), :);
log10isi_nonsomatic = log10isi_all(find(non_somatic), :);

depth_somatic = depth(find(~non_somatic));
depth_nonsomatic = depth(find(non_somatic));

unitType_somatic = unitType(find(~non_somatic));
unitType_nonsomatic = unitType(find(non_somatic));

[WaveformReduced, umap, clusterIdentifiers]=run_umap(waveforms_somatic,...
    'min_dist', 0.1, 'n_neighbors', 20, 'n_components', 2, ...
    'n_epochs', 5000, ...
    'metric', 'euclidean');

[ISIReduced, umap, clusterIdentifiers]=run_umap(log10isi_somatic,...
    'min_dist', 0.1, 'n_neighbors', 20, 'n_components', 2, ...
    'n_epochs', 5000, ...
    'metric', 'euclidean');

%%

NumClusters = 5; % there is a GMM method on Abbaspour's paper; how to find the number of optimal clusters

Attributes = [ISIReduced, WaveformReduced];
CellTypeID = spectralcluster(Attributes, NumClusters, 'Distance', 'mahalanobis'); %'mahalanobis'
% CellTypeID_all = zeros(1, size(waveformsMed_all_norm, 1));
% CellTypeID_all(find(~non_somatic)) = CellTypeID;
% CellTypeID_all(find(non_somatic)) = repmat(NumClusters + 1, 1, length(find(non_somatic)));
% CellTypeID = CellTypeID_all;
% waveforms_somatic = waveformsMed_all_norm;
% log10isi_somatic = log10isi_all;
% depth_somatic = depth;
% unitType_somatic = unitType;

% CellTypeID = unitType_somatic;

% Define unique labels and colormap
unique_labels = unique(CellTypeID);
colors = lines(length(unique_labels)); % Generate distinct colors

% Initialize figure
figure;
% NumClusters = NumClusters + 1;
hold on;
sgtitle(['Number of all cells: ', num2str(length(spikes_somatic))])
% Loop through each label
for i = 1:length(unique_labels)

    subplot(5, NumClusters, i)

    label = unique_labels(i);

    % Extract waveforms corresponding to the current label
    label_waveforms = movmean(waveforms_somatic(CellTypeID == label, :), 1);

    % Compute mean and standard deviation
    avg_waveform = mean(label_waveforms, 1);
    std_waveform = std(label_waveforms, 0, 1);

    % Timepoints for x-axis
    timepoints = linspace(-0.5, 1, 151);

    % Plot the mean waveform

    plot(timepoints, label_waveforms', 'color', [.7 .7 .7])
    hold on
    plot(timepoints, avg_waveform, 'Color', colors(i, :), 'LineWidth', 3, 'DisplayName', sprintf('Label %d', label));

    % Plot shaded region for margins (mean ± std)
    fill([timepoints, fliplr(timepoints)], ...
        [avg_waveform + std_waveform, fliplr(avg_waveform - std_waveform)], ...
        colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylim([0, 1])

    subplot(5, NumClusters, NumClusters+i); hold on

    x_range = -3:0.04:1; % Log10(ISI) range
    celltype = find(CellTypeID == label);
    for celltype_idx = 1:length(celltype)

        plot(10.^x_range, log10isi_somatic(celltype(celltype_idx), :), 'LineWidth', 1.5, 'Color', [.7 .7 .7]);

    end

    plot(10.^x_range, mean(log10isi_somatic(celltype, :)), 'Color', colors(i, :), 'LineWidth', 3, 'DisplayName', sprintf('Label %d', label));

    % Add labels, title, and grid
    xlabel('Log10(ISI)');
    ylabel('Density');
    set(gca, 'XScale', 'log');
    grid on;

    subplot(5, NumClusters, [2*NumClusters+i, 3*NumClusters+i])

    % Scatter plot with specified color
    scatter(rand(1, length(depth_somatic(celltype))) * 5, depth_somatic(celltype), ...
        50, colors(i, :), 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;

    % Box plot adjustments
    box_data = depth_somatic(celltype); % Extract data for the box plot
    box_position = 6; % Position for the boxplot
    box_width = 0.5; % Boxplot width

    % Calculate Q1 and Q3 (for the box limits)
    q1 = prctile(box_data, 25); % First quartile (25th percentile)
    q3 = prctile(box_data, 75); % Third quartile (75th percentile)

    % Define patch vertices for the filled box
    x_patch = [box_position - box_width / 2, box_position + box_width / 2, ...
        box_position + box_width / 2, box_position - box_width / 2];
    y_patch = [q1, q1, q3, q3];

    % Plot the filled box as a patch
    patch(x_patch, y_patch, colors(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Add the actual boxplot lines
    boxchart(ones(size(depth_somatic(celltype))) * 6, depth_somatic(celltype), ...
        'BoxFaceColor', colors(i, :), 'MarkerStyle', 'none', 'LineWidth', 1.5);
    xlim([0, 7])
    title(num2str(length(celltype)))
    ylim([-300, 7500])
    set(gca, 'YDir', 'reverse')
    ylabel('Depth (um)')


    subplot(5, NumClusters, 4*NumClusters+i)

    % Define types corresponding to the numbers
    types = {'NL', 'V', 'S', 'VS', 'R', 'VR', 'SR', 'VSR', 'O'};

    % Count occurrences of each type
    unique_numbers = 0:8; % Numbers 0 to 8
    counts = histcounts(unitType_somatic(celltype), [unique_numbers, max(unique_numbers)+1]); % Histogram for counts

    % Calculate percentages
    total = sum(counts);
    percentages = (counts / total) * 100;

    % Create labels with type names and percentages
    labels = arrayfun(@(idx) sprintf('%s (%.1f%%)', types{idx}, percentages(idx)), ...
        1:length(types), 'UniformOutput', false);

    % Plot pie chart
    pie(counts, labels);

end

hold off;

figure
for i = 1:length(unique_labels)
    label = unique_labels(i);
    mn = mean(movmean(waveforms_somatic(CellTypeID == label, :), 1), 1);
    plot(mn - mn(1));
    hold on
end
