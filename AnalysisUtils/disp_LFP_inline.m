function t = disp_LFP_inline(UIAxes, X, start, ends, offset,feq,ElecName,titre, Col)
% function t = disp_eeg(X,offset,feq,ElecName,titre)
%
% inputs
%     X: dynamics to display. (nbchannels x nbsamples) matrix
%     offset: offset between channels (default max(abs(X)))
%     feq: sapling frequency (default 1)
%     ElecName: cell array of electrode labels (default {S1,S2,...})
%     titre: title of the figure
%
% output
%     t: time vector
%
% G. Birot 2010-02


%% Check arguments
[N K] = size(X);

if nargin < 6
    for n = 1:N
        ElecName{n}  = ['S',num2str(n)];
    end
    titre = [];
end

if nargin < 7
    titre = [];
end

if isempty(feq)
    feq = 1;
end

if isempty(ElecName)
    for n = 1:N
        ElecName{n}  = ['S',num2str(n)];
    end
end

if isempty(offset)
    offset = max(abs(X(:)));
end


%% Build dynamic matrix with offset and time vector
X = X + repmat(offset*(0:-1:-(N-1))',1,K);
t = (-start:ends)/feq;
graduations = offset*(0:-1:-(N-1))';
shiftvec = N:-1:1;
Ysup = max(X(1,:)) + offset;
Yinf = min(X(end,:)) - offset;
% YLabels = cell(N+2) ElecName(shiftvec)

%% Display
% a1 = axes('YAxisLocation','right');
UIAxes.YTickLabel = ElecName(shiftvec);
UIAxes.YTick = graduations(shiftvec);
UIAxes.FontSize = 10;
% a2 = axes('YTickLabel',ElecName(shiftvec),'YTick',graduations(shiftvec),'FontSize',7);
ylim(UIAxes, [Yinf Ysup]);
box(UIAxes, 'on');
grid(UIAxes, 'on')
hold(UIAxes, 'all');

plot(UIAxes, t,X', Col);
% xlabel(UIAxes, 'Time (ms)','FontSize',10);
ylabel(UIAxes, 'Channels','FontSize',10);
xlim(UIAxes, [-start, ends])
title(UIAxes, titre, 'FontSize', 15);
hold(UIAxes, 'off')
