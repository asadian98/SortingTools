function m = plot_CSD_inline(ax, ax2, CSD_matrix,el_pos,dt,scale_plot,max_plot, timeRange, offset, channelNum)
%plot_CSD(CSD_matrix,scale_plot,max_plot,figure_label)
%
%Used to plot the color-plots of the CSD. The colormap_redblackblue goes
%from limit -clim to clim, decided by: clim = max_plot*scale_plot
%
%CSD_matrix: the matrix to plot [A/m^3]
%scale_plot: if one wants to focus the plot this should be from 0 to 1.
%max_plot: the maximum value to plot
%figure_label: text string, e.g. 'a)'


%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html
unit_scale = 1e-3; % A/m^3 -> muA/mm^3
CSD_matrix = CSD_matrix*unit_scale;
%if nargin<3; max_plot=max(abs(CSD_matrix(:))); else; max_plot = max_plot*unit_scale; end;
%if nargin<2; scale_plot = 1; end;
if max_plot==0; max_plot=max(abs(CSD_matrix(:))); end;

%figure();
clim=max_plot*scale_plot;

npoints = 200; % number of points to plot in the vertical direction
le = length(el_pos);
first_z = el_pos(1)-(el_pos(2)-el_pos(1))/2; %plot starts at z1-h/2;
last_z = el_pos(le)+(el_pos(le)-el_pos(le-1))/2; %ends at zN+h/2;
zs = first_z:(last_z-first_z)/npoints:last_z;
el_pos(le+1) = el_pos(le)+(el_pos(le)-el_pos(le-1)); % need this in for loop
j=1; %counter
for i=1:length(zs) % all new positions
    if zs(i)>(el_pos(j)+(el_pos(j+1)-el_pos(j))/2) % > el_pos(j) + h/2
        j = min(j+1,le);
    end;
    new_CSD_matrix(i,:)=CSD_matrix(j,:);
end;
CSD_matrix = new_CSD_matrix;
%
% npoints = 200; % number of points to plot in the vertical direction
% le = length(el_pos);
% first_z = el_pos(1)-(el_pos(2)-el_pos(1))/2; %plot starts at z1-h/2;
% last_z = el_pos(le)+(el_pos(le)-el_pos(le-1))/2; %ends at zN+h/2;
% zs = first_z:(last_z-first_z)/npoints:last_z;
% el_pos(le+1) = el_pos(le)+(el_pos(le)-el_pos(le-1)); % need this in for loop
% j=1; %counter
% for i=1:length(zs) % all new positions
%     if zs(i)>(el_pos(j)+(el_pos(j+1)-el_pos(j))/2) % > el_pos(j) + h/2
%         j = min(j+1,le);
%     end;
%     new_CSD_matrix(i,:)=CSD_matrix(j,:);
% end;
% % %increase z-resolution:
% % inc_factor = 10; % 10 times the resolution
% % for i = 1:length(CSD_matrix(:,1))*inc_factor
% %     new_CSD_matrix(i,:) = CSD_matrix(round((i-1)/inc_factor+0.5),:);
% % end;
[m, n] = size(CSD_matrix);

imagesc(ax, CSD_matrix,[-clim clim]);
colormap(ax, colormap_redblackblue), colorbar(ax);
set(ax,'FontSize',10)

[nel,ntime] = size(CSD_matrix);
time = dt:dt:dt*ntime;

xticks(ax, 100:100:(-offset(1)+offset(end)-1))

time_ticks = get(ax, 'XTick');
for i=1:length(time_ticks)
    new_ticks{i} = num2str(time(time_ticks(i)) + offset(1));
end

yticks(ax, m/(channelNum*2):m/channelNum:m);

z_ticks = get(ax,'YTick');
for i=channelNum:-1:1
    new_zticks{i} = num2str(i); % [m] -> [mm] with two decimals
end

% xticks(ax, -500:100:1000)
offset(1)
offset(end)
xlim(ax, [0, -offset(1) + offset(end)])
set(ax, 'XTickLabel',new_ticks)
ylim(ax, [0, m])
% set(ax, 'YTick',z_ticks)
set(ax, 'YTickLabel',new_zticks)

% ------------------------------------ AX 2 ------------------------------
CSDprofile = mean(CSD_matrix(:, timeRange), 2);
plot(ax2, CSDprofile, m:-1:1)
yticks(ax2, m/(channelNum*2):m/channelNum:m);
ylim(ax2, [1, m])
z_ticks = get(ax2,'YTick');
tmp = channelNum:-1:1;
for i=1:channelNum
    new_zticks{i} = num2str(tmp(i)); % [m] -> [mm] with two decimals
end
set(ax2, 'YTickLabel',new_zticks, 'box', 'off', 'TickDir', 'out', 'FontSize', 10)
xlabel(ax2, 'CSD (A m^{-2})')
xline(ax2, 0, '--r')

[~, I] = max(-CSDprofile);

tmp = find(CSDprofile(1:I) > 0);
if(~isempty(tmp))
    yline(ax2, m-tmp(end), '--b')
    ch_up = tmp(end)*channelNum/200;
    ch_down = (tmp(end)-1)*channelNum/200+1;
    ch = round((ch_up + ch_down) / 2);
    text(ax2, min(CSDprofile)+5, m-tmp(end)+10, ['Channel ', num2str(ch)]);
else
end

