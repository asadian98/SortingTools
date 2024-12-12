function runKilo(Folders, meta, thresh)

        nCh = meta.nCh;
        S_probe16 = meta.S_probe16;
        Fs = meta.fs;
        sesName = meta.sesName;
        fileName = meta.fileName;
        location = meta.location;

        connected       = true(nCh, 1);
        chanMap         = 1:nCh;
        chanMap0ind     = chanMap - 1;

        % define the channel map as a filename (string) or simply an array
        if(S_probe16)
            probeGeometry  = 'linear300'; % 300 um spacing, change this for Grover
        else
            probeGeometry  = 'linear150'; % 150 um spacing, change this for Grover
        end

        ops.probeGeometry  = probeGeometry;  % see probeGeometry2coords.m for options

        [xcoords, ycoords, kcoords] = probeGeometry2coords(probeGeometry, nCh);

        % save:
        save(fullfile(Folders.sortingFolder, 'chanMap150.mat'), ...
            'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'Fs', 'probeGeometry')

        %%

        ops.root                = [Folders.KiloFolder, sesName];
        ops.chanMap             = fullfile(Folders.sortingFolder, 'chanMap150.mat');
        ops.fshigh = 150;
        ops.minfr_goodchannels = 0.05; % Ks default is 0.1
        ops.Th = [6 4];
        ops.lam = 10;
        % splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
        ops.AUCsplit = 0.92; % 0.92
        % minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
        ops.minFR = 1/50;
        % spatial constant in um for computing residual variance of spike
        ops.sigmaMask = 30;
        % threshold crossings for pre-clustering (in PCA projection space)
        ops.ThPre = 8;

        ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
        ops.parfor              = 0; % whether to use parfor to accelerate some parts of the algorithm
        ops.verbose             = 1; % whether to print command line progress
        ops.showfigures         = 0; % whether to plot figures during optimization

        ops.datatype            = 'dat';  % binary ('dat', 'bin') or 'openEphys'

        ops.fs                  = Fs;        % sampling rate		(omit if already in chanMap file)
        ops.NchanTOT            = nCh;           % total number of channels (omit if already in chanMap file)
        ops.Nchan               = nCh;           % number of active channels (omit if already in chanMap file)
        % ops.Nfilt               = ceil(nCh / 32) * 32 * 4;           % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)

        % options for channel whitening
        ops.whitening           = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)

        % danger, changing these settings can lead to fatal errors
        ops.spkTh           = thresh;      % spike threshold in standard deviations (-6)
        ops.reorder         = 1;       % whether to reorder batches for drift correction.
        ops.nskip           = 25;  % how many batches to skip for determining spike PCs

        % ops.Nfilt               = 1024; % max number of clusters
        ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
        ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
        ops.NT                  = 64*1024 + ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory).
        ops.whiteningRange      = 32; % number of channels to use for whitening each channel (24?)
        ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch (25?)
        ops.scaleproc           = 200;   % int16 scaling of whitened data
        ops.nPCs                = 3; % how many PCs to project the spikes into
        ops.useRAM              = 0; % not yet available
        % number of samples to average over (annealed from first to second value)
        ops.momentum = [20 400];

        %% this block runs all the steps of the algorithm
        rootZ = [Folders.KiloFolder, sesName]; % the raw data binary file is in this folder
        rootH = [Folders.KiloFolder, 'temp\']; % path to temporary binary file (same size as data, should be on fast SSD)

        ops.trange      = [0 Inf]; % time range to sort

        ops.fproc       = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD

        fprintf('Looking for data inside %s \n', rootZ)

        % is there a channel map file in this folder?
        fs = dir(fullfile(rootZ, 'chan*.mat'));
        if ~isempty(fs)
            ops.chanMap = fullfile(rootZ, fs(1).name);
        end

        % find the binary file
        fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
        ops.fbinary = fullfile(rootZ, fs(1).name);

        % preprocess data to create temp_wh.dat
        rez = preprocessDataSub(ops);

        % time-reordering as a function of drift
        rez = clusterSingleBatches(rez);

        % % saving here is a good idea, because the rest can be resumed after loading rez
        % save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

        % main tracking and template matching algorithm
        rez = learnAndSolve8b(rez);

        % final merges
        rez = find_merges(rez, 1);

        % final splits by SVD
        rez = splitAllClusters(rez, 1);

        % final splits by amplitudes
        rez = splitAllClusters(rez, 0);

        % decide on cutoff
        rez = set_cutoff(rez);

        close all

        fprintf('found %d good units \n', sum(rez.good>0))

        % write to Phy
        fprintf('Saving results to Phy  \n')
        rezToPhy(rez, Folders.sortingFolder);
        %
        % % if you want to save the results to a Matlab file...
        %
        % % discard features in final rez file (too slow to save)
        % rez.cProj = [];
        % rez.cProjPC = [];
        %
        % % save final results as rez2
        % fprintf('Saving final results in rez2  \n')
        % fname = fullfile(rootZ, 'rez2.mat');
        % save(fname, 'rez', '-v7.3');

        % meta info:
        info.dsn            = sesName;
        info.rawFolder      = Folders.Ripple_Raw_Data_Folder;
        info.rawFile        = fileName;
        info.datestr        = datestr(now, 'yyyymmddTHHMM');
        info.RecLocation    = location;
        info.ops = ops;
        save(fullfile(Folders.sortingFolder, 'convertInfo.mat'),  'info');
        disp('info - saved')

        % sampsToSecsMap:
        % transform spike times in samples to spike times in seconds and store.

        cd(Folders.sortingFolder)
        spikeTimesSamples = readNPY(fullfile(Folders.sortingFolder, ...
            'spike_times.npy'));
        spikeTimesSeconds = spikeTimesSamples/Fs;
        disp('Saving spike times in seconds to .npy file.')
        writeNPY(spikeTimesSeconds, fullfile(Folders.sortingFolder, ...
            'spike_times_seconds.npy'));

        disp('Enjoy working with the data! Phy will pop up shortly ... ')

end