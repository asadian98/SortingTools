function my_csdViewer(cont_tr, st, en, varargin)
% function my_csdViewer(spikeTimes, clu, eventTimes, window, trGroups[, params])
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
fprintf(1, '- left/right arrow: select previous/next Trial\n')
fprintf(1, '- up/down arrow: change smoothing of csd\n')
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


tr = 1:size(cont_tr, 2);
CSD_cs = {};
for j = 1:length(tr)

    lfp_data = reshape(cont_tr(:, tr(j), :), size(cont_tr(:, 1, :), 1), []);
    channelNum = size(lfp_data, 1);

    % filter parameters:
    gauss_sigma = 0.001*1e-3; % mm -> m
    filter_range = 5*gauss_sigma; % numeric filter must be finite in extent
    % electrical parameters:
    cond = 0.3;
    cond_top = 0.3;

    % size, potential (m1 has to equal number of electrode contacts)
    m1 = size(lfp_data, 1);

    % geometrical parameters:
    diam = 0.5*1e-3;

    ch_num = size(cont_tr, 1);
    if(ch_num == 32)
        step = 0.15;
    else
        step = 0.3;
    end
    el_pos = [0.1:step:4.75]*1e-3;

    % compute spline iCSD:
    Fcs = F_cubic_spline(el_pos,diam,cond,cond_top);
    [zs,CSD_cs{j}] = make_cubic_splines(el_pos, lfp_data(:, :), Fcs);

    if gauss_sigma~=0 %filter iCSD
        [zs,CSD_cs{j}]=gaussian_filtering(zs,CSD_cs{j},gauss_sigma,filter_range);
    end

    b1 = 0.2300;

    b0 = 0.5400;

    % compute delta iCSD:
    CSD_cs{j} = F_delta(el_pos,diam,cond,cond_top)^-1*lfp_data(:, :);
    
    if b1~=0 %filter iCSD
      [n1,n2]=size(CSD_cs{j});            
      CSD_add(1,:) = zeros(1,n2);   %add top and buttom row with zeros
      CSD_add(n1+2,:)=zeros(1,n2);
      CSD_add(2:n1+1,:)=CSD_cs{j};        %CSD_add has n1+2 rows
      CSD_cs{j} = S_general(n1+2,b0,b1)*CSD_add; % CSD has n1 rows
    end

    zs = el_pos;

end

f = figure; f.Color = 'w';
set(f, 'KeyPressFcn', @(f,k)csdViewerCallback(f, k));

subplot(2, 2, 1)
ch_Arr = [];
CSD = [];
for j = 1:length(CSD_cs)

    CSDprofile = mean(CSD_cs{j}(:, 500+(st:en)), 2);
    [m, n] = size(CSD_cs{j});

    plot(CSDprofile, m:-1:1, 'k')
    hold on;
    z_ticks = get(gca,'YTick');
    yticks = linspace(z_ticks(1), z_ticks(end), channelNum);
    yticklabels = channelNum:-1:1;
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

    CSD = [CSD, CSDprofile];

    [~, I] = max(-CSDprofile);
    tmp = find(CSDprofile(1:I) > 0);
    if(~isempty(tmp))
        ch_up = tmp(end)*channelNum/200;
        ch_down = (tmp(end)-1)*channelNum/200+1;
        ch = round((ch_up + ch_down) / 2);
        ch_Arr = [ch_Arr, ch];
    end
    
end

disp(['The most frequent number is: ', num2str(mode(ch_Arr))]);

plot(mean(CSD, 2), m:-1:1, 'r', 'LineWidth',2)
CSDprofile = mean(CSD, 2);
[~, I] = max(-CSDprofile);

tmp = find(CSDprofile(1:I) > 0);
if(~isempty(tmp))
    yline(gca, m-tmp(end), '--b')
    ch_up = tmp(end)*channelNum/200;
    ch_down = (tmp(end)-1)*channelNum/200+1;
    ch = round((ch_up + ch_down) / 2);
    title(['Channel ', num2str(ch)]);
end

myData.CSD_cs = CSD_cs;
myData.cont_tr = cont_tr;
myData.tr = tr;
myData.params.tr_want = 1;
myData.trialNum = size(cont_tr, 2);
myData.zs = zs;
myData.channelNum = channelNum;
myData.st = st;
myData.en = en;

set(f, 'UserData', myData);

% Make plots
myData.plotAxes = [];
if isempty(myData.plotAxes)
    for pidx = 1:4
        subplot(2,2,pidx);
        myData.plotAxes(pidx) = gca;
    end
    set(f, 'UserData', myData);
end

csdViewerPlot(f)

end

function csdViewerPlot(f)

myData = get(f,'UserData');
cont_tr = myData.cont_tr;
CSD_cs = myData.CSD_cs;
tr = myData.tr;
zs = myData.zs;
p = myData.params;
tr_want = p.tr_want;
channelNum = myData.channelNum;
st = myData.st;
en = myData.en;

subplot(2, 2, 2)
cla(gca)
disp_LFP_inline(gca, reshape(cont_tr(:, tr(tr_want), :), size(cont_tr(:, 1, :), 1), []), 500, 500, 100, 1, [], '', 'k');
title(['Trial: ', num2str(tr_want)])
x = [myData.st myData.st myData.en myData.en];
yl = ylim;
y = [yl(1) yl(2) yl(2) yl(1)];
patch(x,y,'red','FaceAlpha',.3);


subplot(2, 2, 4)
cla(gca)
dt = 0.5;
scale_plot = 1;
max_plot = 0;
plot_CSD(CSD_cs{tr_want},zs,dt,scale_plot,max_plot)
colorbar off
xlim([0, 1000])
x_ticks = get(gca,'XTick');
xticks = linspace(0, 1000, 5);
xticklabels = linspace(-500, 500, 5);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
colorbar

subplot(2, 2, 3)
cla(gca)
[m, n] = size(CSD_cs{tr_want});
CSDprofile = mean(CSD_cs{tr_want}(:, 500+(st:en)), 2);
plot(CSDprofile, m:-1:1, 'k')
hold on;
z_ticks = get(gca,'YTick');
yticks = linspace(z_ticks(1), z_ticks(end), channelNum);
yticklabels = channelNum:-1:1;
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

[~, I] = max(-CSDprofile);

tmp = find(CSDprofile(1:I) > 0);
if(~isempty(tmp))
    yline(gca, m-tmp(end), '--b')
    ch_up = tmp(end)*channelNum/200;
    ch_down = (tmp(end)-1)*channelNum/200+1;
    ch = round((ch_up + ch_down) / 2);
    text(gca, min(CSDprofile)-1, m-tmp(end)+10, ['Channel ', num2str(ch)]);
end
xlabel('CSD (A m^{-2})')
xline(0)

end

function csdViewerCallback(f, keydata)

% fprintf('callback on %d with source %d\n', f.Number, keydata.Source.Number);

myData = get(f, 'UserData');
p = myData.params;

switch keydata.Key
    case 'rightarrow' % increment cluster index

        p.tr_want = p.tr_want+1;
        if p.tr_want>myData.trialNum
            p.tr_want=1;
        end
       
    case 'leftarrow' % decrement cluster index

        p.tr_want = p.tr_want-1;
        if p.tr_want<1
            p.tr_want=myData.trialNum;
        end

    case 'uparrow' % increase smoothing
%         p.smoothSize = p.smoothSize*1.2;
%         p.smWin = genSmWin(p);
        myData.st = myData.st + 50;
        myData.en = myData.en + 50;
        
    case 'downarrow' % decrease smoothing
%         p.smoothSize = p.smoothSize/1.2;
%         p.smWin = genSmWin(p);

        myData.st = myData.st - 50;
        myData.en = myData.en - 50;

    case 'f'
        p.filterType = p.filterType+1;
        if p.filterType>3; p.filterType = 1; end

        p.smWin = genSmWin(p);

    case 'e' % whether to show standard error as shading
        p.showErrorShading = ~p.showErrorShading;

    case 't' % whether to plot the psth trace for each condition or just the overall one
        p.showAllTraces = ~p.showAllTraces;



    case 'r'
        ax = subplot(3,1,1); title('click start and stop of range')
        %         [startRange,~] = ginput();
        %         [stopRange,~] = ginput();
        waitforbuttonpress;
        q = get(ax, 'CurrentPoint');
        p.startRange = q(1,1);
        waitforbuttonpress;
        q = get(ax, 'CurrentPoint');
        p.stopRange = q(1,1);
        if p.stopRange<p.startRange
            tmp = p.startRange;
            p.startRange = p.stopRange;
            p.stopRange = tmp;
        end

    case 'c'
        newC = inputdlg('cluster ID?');
        ind = find(myData.clusterIDs==str2num(newC{1}),1);
        if ~isempty(ind)
            p.clusterIndex = ind;
        end

    case 'x'
        p.colorType = p.colorType+1;
        p.colors = genColors(p.colorType, myData.nGroups);

    case '2' % increase tick size
        p.rasterScale = p.rasterScale*1.5;

    case '1' % decrease tick size
        p.rasterScale = p.rasterScale/1.5;
end

myData.params = p;
set(f, 'UserData', myData);

% plot with new settings
csdViewerPlot(f)


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
function exitCallback(src, event, f)
myData = get(f, 'UserData');
su = myData.su;
Z_Name = myData.params.Z_Name;
myData.Z_Spikes.su = su;
Z_Spikes = myData.Z_Spikes;
save(Z_Name, 'Z_Spikes')
disp('Saved')
end