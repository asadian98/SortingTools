function plot_rasters_above(SPIKES,color,yaxis_add,optionalmarks1,optionalmarks2)
% Simple function to add rasters above spikedensity functions
% Inputs:
%   Spikes: Large arrays of NaNs and spike times, relatively to current xaxis
%   color: Desired tic colour
%   yaxis_add: How much above yaxis to add
%   optionalmarks1 and optionalmarks2: self explanatory
XAXIS = xlim;
YAXIS = ylim;   
YAXIS_start = YAXIS(2);
YINC = yaxis_add/size(SPIKES,1);
TIC_HEIGHT = YINC * .9;
for a = 1:size(SPIKES,1)    % For every trial
   spikes = SPIKES(a,:);
   spikes = round(spikes(~isnan(spikes)));
   YAXIS_start = YAXIS_start + YINC;
   for x = 1:length(spikes)
       plot([spikes(x) spikes(x)],[YAXIS_start YAXIS_start+TIC_HEIGHT],'LineWidth',2,'Color',color);
   end
   if nargin > 3
       plot(optionalmarks1(a),YAXIS_start+(YINC/2),'ks','MarkerFaceColor','k', 'LineWidth',0.1);
   end
   if nargin > 4
       plot(optionalmarks2(a),YAXIS_start+(YINC/2), 'o','MarkerEdgeColor', [0.9290 0.6940 0.1250],'MarkerFaceColor', [0.9290 0.6940 0.1250]	, 'LineWidth',0.5);
   end

end
axis([XAXIS(1) XAXIS(2) 0 YAXIS_start+YINC]);

