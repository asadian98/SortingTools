function my_eyetraceViewer(Z_ML, varargin)
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
fprintf(1, '- r: change to rightward saccade\n')
fprintf(1, '- l: change to leftward saccade\n')
fprintf(1, '- u: change to upward saccade\n')
fprintf(1, '- d: change to downward saccade\n')
fprintf(1, '- f: saccade is fix in RL\n')
fprintf(1, '- g: saccade is fix in UD\n')
fprintf(1, '- s: save Z_ML\n')

fprintf(1, '- 1/2: decrease/increase the raster tick size\n')

if ~isempty(varargin)
    p = varargin{1}; % parameters supplied by user
else
    p = [];
end

myData.params.Z_Name = getOr(p, 'Z_Name', 'HI');

f = figure; f.Color = 'w';
set(f, 'KeyPressFcn', @(f,k)eyetraceCallback(f, k));

if(~isfield(Z_ML, 'EyeDir_All_RL'))
Z_ML.EyeRT_All = [];
Z_ML.EyeDir_All_RL = [];
Z_ML.EyeDir_All_UD = [];

for i = 1:length(Z_ML.condition)

    if Z_ML.TrialError(i) > 0 % Only run saccade detection on "correct" trials
        Eh = Z_ML.EyeX_Filt_Diode{i}; % Horizontal Eye Position
        Ev = Z_ML.EyeY_Filt_Diode{i}; % Vertical Eye Position
        dEh = diff(Eh')*1000;
        dEv = diff(Ev')*1000;

        DiodeTime = 0; % I want to detect all saccades even before the PD
        output=sacdetector(dEh,dEv,Eh',Ev',[100 50 50],DiodeTime);
        if ~isempty(output)
            Z.movearray{i}=output;

            AAA = output(:,3);   % RTs relative to target onset
            AAA = AAA(AAA>0); % RTs for eye movements occurring after target onset
            RT = AAA;
            RT = RT(find(RT < length(Z_ML.EyeX_Filt_Diode{i})-50));
            RT = RT(find(RT > 11));

            RL = zeros(1, length(RT));
            UD = zeros(1, length(RT));
            for j = 1:length(RT)
                x_diff = Z_ML.EyeX_Filt_Diode{i}(DiodeTime + RT(j) + 50) - Z_ML.EyeX_Filt_Diode{i}(DiodeTime + RT(j) -10);
                y_diff = Z_ML.EyeY_Filt_Diode{i}(DiodeTime + RT(j) + 50) - Z_ML.EyeY_Filt_Diode{i}(DiodeTime + RT(j) -10);

                if(abs(x_diff)<2)  RL(j) = 3; elseif(x_diff > 0) RL(j) = 1; elseif(x_diff < 0) RL(j) = 2; end
                if(abs(y_diff)<2)  UD(j) = 3; elseif(y_diff > 0) UD(j) = 1; elseif(y_diff < 0) UD(j) = 2; end
            end

            idxx = find(RL == 3 & UD == 3);
            RL(idxx) = [];
            UD(idxx) = [];
            RT(idxx) = [];

            Z_ML.EyeRT_All(i, 1:length(RT)) = RT;
            Z_ML.EyeDir_All_RL(i, 1:length(RL)) = RL;
            Z_ML.EyeDir_All_UD(i, 1:length(UD)) = UD;
        end
    end
end
end
% Create a text box
hTextBox = uicontrol('Style', 'edit', 'Position', [1000 20 80 40], 'String', 'RL');
% Create a text box
hTextBox2 = uicontrol('Style', 'edit', 'Position', [1200 20 80 40], 'String', 'UP');

% Add text box data to the figure's UserData
myData.hTextBox = hTextBox;
myData.textBoxData = '';
set(f, 'UserData', myData);

% Add text box data to the figure's UserData
myData.hTextBox2 = hTextBox2;
myData.textBoxData2 = '';
set(f, 'UserData', myData);

myData.Z_ML = Z_ML;
myData.params.tr = find(Z_ML.TrialError == 1, 1);
myData.params.sac_idx = 1;
set(f, 'UserData', myData);

eyetraceViewerPlot(f)

end

function eyetraceViewerPlot(f)

cla
myData = get(f,'UserData');
Z_ML = myData.Z_ML;
tr = myData.params.tr;
sac_idx = myData.params.sac_idx;

subplot(1, 2, 1)
plot(Z_ML.EyeX_Filt_Diode{tr}, 'b')
xline(Z_ML.EyeRT_All(tr, :), 'k')
xline(Z_ML.EyeRT_All(tr, sac_idx), 'r', 'LineWidth', 3)
view([90 -90]);

sgtitle(['Trial number:     ', num2str(tr), '     b: X,    r: Y'])

subplot(1, 2, 2)
plot(Z_ML.EyeY_Filt_Diode{tr}, 'r')
hold on
xline(Z_ML.EyeRT_All(tr, :), 'k')
xline(Z_ML.EyeRT_All(tr, sac_idx), 'r', 'LineWidth', 3)



if(Z_ML.EyeDir_All_RL(tr, sac_idx) == 1)
    set(myData.hTextBox, 'String', 'Right');
elseif(Z_ML.EyeDir_All_RL(tr, sac_idx) == 2)
    set(myData.hTextBox, 'String', 'Left');
elseif(Z_ML.EyeDir_All_RL(tr, sac_idx) == 3)
    set(myData.hTextBox, 'String', 'Fix');
else
    set(myData.hTextBox, 'String', 'no saccade');
end

if(Z_ML.EyeDir_All_UD(tr, sac_idx) == 1)
    set(myData.hTextBox2, 'String', 'UP');
elseif(Z_ML.EyeDir_All_UD(tr, sac_idx) == 2)
    set(myData.hTextBox2, 'String', 'Down');
elseif(Z_ML.EyeDir_All_UD(tr, sac_idx) == 3)
    set(myData.hTextBox2, 'String', 'Fix');
else
    set(myData.hTextBox2, 'String', 'no saccade');
end


end

function eyetraceCallback(f, keydata)

% fprintf('callback on %d with source %d\n', f.Number, keydata.Source.Number);

myData = get(f, 'UserData');
Z_ML = myData.Z_ML;
p = myData.params;

switch keydata.Key
    case 'rightarrow' % increment cluster index

        idx = find(Z_ML.TrialError == 1);
        if(find(idx == p.tr) ~= length(idx))
            p.tr = idx(find(idx == p.tr) + 1);
        else
            p.tr= idx(1);
        end

        p.sac_idx = 1;

    case 'leftarrow' % decrement cluster index

        idx = find(Z_ML.TrialError == 1);
        if(find(idx == p.tr) ~= 1)
            p.tr = idx(find(idx == p.tr) - 1);
        else
            p.tr= idx(end);
        end
        p.sac_idx = 1;
    case 'uparrow' % go to next saccade

        idx = find(Z_ML.EyeRT_All(p.tr, :) ~= 0);
        if(find(idx == p.sac_idx) ~= length(idx))
            p.sac_idx = idx(find(idx == p.sac_idx) + 1);
        else
            p.sac_idx = idx(1);
        end

    case 'downarrow' % decrease smoothing
        idx = find(Z_ML.EyeRT_All(p.tr, :) ~= 0);
        if(find(idx == p.sac_idx) ~= 1)
            p.sac_idx = idx(find(idx == p.sac_idx) - 1);
        else
            p.sac_idx= idx(end);
        end

        %     case 'f'
        %         %         p.filterType = p.filterType+1;
        %         %         if p.filterType>3; p.filterType = 1; end
        %         %
        %         %         p.smWin = genSmWin(p);
        %
        %     case 'e' % whether to show standard error as shading
        %         %         p.showErrorShading = ~p.showErrorShading;
        %
        %     case 't' % whether to plot the psth trace for each condition or just the overall one
        %         p.showAllTraces = ~p.showAllTraces;

    case 'b'
        i = p.tr;
        if Z_ML.TrialError(i) > 0 % Only run saccade detection on "correct" trials
            Eh = Z_ML.EyeX_Filt_Diode{i}; % Horizontal Eye Position
            Ev = Z_ML.EyeY_Filt_Diode{i}; % Vertical Eye Position
            dEh = diff(Eh')*1000;
            dEv = diff(Ev')*1000;

            DiodeTime = 0; % I want to detect all saccades even before the PD
            output=sacdetector(dEh,dEv,Eh',Ev',[100 50 50],DiodeTime);
            if ~isempty(output)
                Z.movearray{i}=output;

                AAA = output(:,3);   % RTs relative to target onset
                AAA = AAA(AAA>0); % RTs for eye movements occurring after target onset
                RT = AAA;
                RT = RT(find(RT < length(Z_ML.EyeX_Filt_Diode{i})-50));
                RT = RT(find(RT > 11));

                RL = zeros(1, length(RT));
                UD = zeros(1, length(RT));
                for j = 1:length(RT)
                    x_diff = Z_ML.EyeX_Filt_Diode{i}(DiodeTime + RT(j) + 50) - Z_ML.EyeX_Filt_Diode{i}(DiodeTime + RT(j) -10);
                    y_diff = Z_ML.EyeY_Filt_Diode{i}(DiodeTime + RT(j) + 50) - Z_ML.EyeY_Filt_Diode{i}(DiodeTime + RT(j) -10);

                    if(abs(x_diff)<2)  RL(j) = 3; elseif(x_diff > 0) RL(j) = 1; elseif(x_diff < 0) RL(j) = 2; end
                    if(abs(y_diff)<2)  UD(j) = 3; elseif(y_diff > 0) UD(j) = 1; elseif(y_diff < 0) UD(j) = 2; end
                end

                idxx = find(RL == 3 & UD == 3);
                RL(idxx) = [];
                UD(idxx) = [];
                RT(idxx) = [];

                Z_ML.EyeRT_All(i, 1:length(RT)) = RT;
                Z_ML.EyeDir_All_RL(i, 1:length(RL)) = RL;
                Z_ML.EyeDir_All_UD(i, 1:length(UD)) = UD;
            end
        end

    case 'c' % Remove a saccade

        if(p.sac_idx < length(Z_ML.EyeRT_All(p.tr, :)))
            Z_ML.EyeRT_All(p.tr, p.sac_idx:length(Z_ML.EyeRT_All(p.tr, :))) = [Z_ML.EyeRT_All(p.tr, (p.sac_idx+1):length(Z_ML.EyeRT_All(p.tr, :))), 0];
            Z_ML.EyeDir_All_RL(p.tr, p.sac_idx:length(Z_ML.EyeDir_All_RL(p.tr, :))) = [Z_ML.EyeDir_All_RL(p.tr, (p.sac_idx+1):length(Z_ML.EyeDir_All_RL(p.tr, :))), 0];
            Z_ML.EyeDir_All_UD(p.tr, p.sac_idx:length(Z_ML.EyeDir_All_UD(p.tr, :))) = [Z_ML.EyeDir_All_UD(p.tr, (p.sac_idx+1):length(Z_ML.EyeDir_All_UD(p.tr, :))), 0];
        else
            Z_ML.EyeRT_All(p.tr, p.sac_idx) = 0;
            Z_ML.EyeDir_All_RL(p.tr, p.sac_idx) = 0;
            Z_ML.EyeDir_All_UD(p.tr, p.sac_idx) = 0;
        end

        idx = find(Z_ML.EyeRT_All(p.tr, :) ~= 0);
        if(find(idx == p.sac_idx) < length(idx))

        else
            p.sac_idx = idx(1);
        end

    case 'x'
        newC = inputdlg('trial number?');
        p.tr = str2num(newC{1});

    case 'r'
        Z_ML.EyeDir_All_RL(p.tr, p.sac_idx) = 1;

    case 'l'
        Z_ML.EyeDir_All_RL(p.tr, p.sac_idx) = 2;
        
    case 'u'
        Z_ML.EyeDir_All_UD(p.tr, p.sac_idx) = 1;

    case 'd'
        Z_ML.EyeDir_All_UD(p.tr, p.sac_idx) = 2;

    case 'f'
        Z_ML.EyeDir_All_RL(p.tr, p.sac_idx) = 3;

    case 'g'
        Z_ML.EyeDir_All_UD(p.tr, p.sac_idx) = 3;

    case 's'
        save(p.Z_Name, 'Z_ML');
        msgbox('Saved!', 'Saved');

end

myData.params = p;
myData.Z_ML = Z_ML;
set(f, 'UserData', myData);

% plot with new settings
eyetraceViewerPlot(f)


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