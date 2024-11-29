function plot_rasters_above_app(fig_handle, SPIKES,color,yaxis_add,optionalmarks1,optionalmarks2)
% Simple function to add rasters above spikedensity functions
% Inputs:
%   Spikes: Large arrays of NaNs and spike times, relatively to current xaxis
%   color: Desired tic colour
%   yaxis_add: How much above yaxis to add
%   optionalmarks1 and optionalmarks2: self explanatory
XAXIS = xlim(fig_handle);
YAXIS = ylim(fig_handle);   
YAXIS_start = YAXIS(2);
YINC = yaxis_add/size(SPIKES,1);
TIC_HEIGHT = YINC * .9;
for a = 1:size(SPIKES,1)    % For every trial
   spikes = SPIKES(a,:);
   spikes = round(spikes(~isnan(spikes)));
   YAXIS_start = YAXIS_start + YINC;
   for x = 1:length(spikes)
       plot(fig_handle, [spikes(x) spikes(x)],[YAXIS_start YAXIS_start+TIC_HEIGHT],'LineWidth',1.5,'Color',color);
   end
   if nargin > 4
       plot(fig_handle, optionalmarks1(a),YAXIS_start+(YINC/2),'ks','MarkerFaceColor','k');
   end
   if nargin > 5
       plot(fig_handle, optionalmarks2(a),YAXIS_start+(YINC/2),'ro','MarkerFaceColor','r');
   end

end
axis(fig_handle, [XAXIS(1) XAXIS(2) 0 YAXIS_start+YINC]);

