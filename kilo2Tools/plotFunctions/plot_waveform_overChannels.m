function hA = plot_waveform_overChannels(sp, flag, cid)
hold off
chList          = 1:sp.n_channels_dat;
chListFlip      = fliplr(chList);
chListFlipNeg   = -chListFlip;
yScaler         = 0.01;

for iS = find(sp.clusterScore == 2 | sp.clusterScore == 1) % Just plot 'good' + 'MU' clusters
    ch = double(sp.peakCh(iS));
    wf = sp.medWfOnPeakCh(iS,:) .* yScaler;
    wfLength = numel(wf);

    x = linspace(-0.5, 1.25, wfLength);
%     x  = (1:wfLength)  - ceil(wfLength/2);
    if(sp.clusterId(iS) == cid-1)
        hP(iS) = plot(x, -ch + wf, 'r', 'LineWidth', 3); % flipping sign for plotting
        title(['Waveform -- Score: ', num2str(sp.clusterScore(iS))])
    else
        if(sp.clusterScore(iS) == 2)
            hP(iS) = plot(x, -ch + wf, 'g'); % flipping sign for plotting
        else
            hP(iS) = plot(x, -ch + wf, 'k'); % flipping sign for plotting
        end
    end
    hold on
    xJitter = 0;
    yJitter = 0.1;
    if(~flag)
        text(-0.5/2 + xJitter, -ch + yJitter, ['cid ' num2str(sp.clusterId(iS)+1)], 'Color', 'k', 'FontSize', 6)
    end
end 

ylim([chListFlipNeg(1)-1 chListFlipNeg(end)+1])
xlim(1.25.*[-0.5 1.25])
ylabel('Channel Number')
xlabel('time(ms)')
set(gca, 'YTick', chListFlipNeg, 'YTickLabel', chListFlip)
