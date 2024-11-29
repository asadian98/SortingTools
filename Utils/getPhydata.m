function getPhydata(ses_name, datFolder, region)

cd(datFolder)

% Maybe useful: https://phy.readthedocs.io/en/latest/sorting_user_guide/

spike_times = readNPY('spike_times.npy');
spike_clusters = readNPY('spike_clusters.npy');
% load('sp')
% spike_templates = readNPY('spike_templates.npy');

% Initially, before running phy, the spike-cluster and spike-template
% assignments are identical. If spike_clusters.npy does not exist,it is 
% automatically copied from spike_templates.npy. When modifying the spike-cluster assignments in phy, 
% only spike_clusters.npy is modified, while spike_templates.npy remains
% unchanged. --> use spike_clusters

% t = readtable("cluster_group.tsv", "FileType","delimitedtext");
t2 = readtable("cluster_channel.tsv", "FileType","delimitedtext");
chanMap = readNPY('channel_map.npy');
chanMap = chanMap + 1;
% cluster_group = [];
% 
% for i = 1:size(t, 1)
%     
%    cluster_group(i, 1) = t{i, 1} + 1;
%    if(strcmp(t{i, 2}, 'good'))
%        cluster_group(i, 2) = 1;
%    else
%        cluster_group(i, 2) = 0;
%    end
%     
% end
% 
% cluster_group = cluster_group';
spikes_phy = [double(spike_clusters)+1, double(spike_times)]'; 
csvFile         = [datFolder, '\cluster_group.tsv'];
[clusterId, clusterScore]  = readClusterGroupsCSV(csvFile);

cluster_channel = table2array(t2(:, 2)) + 1;

% Remove noise clusters
% idx_notNoise = find(clusterScore ~= 0);
clusterId = clusterId + 1;
cluster_channel = cluster_channel;
% cluster_channel = sp.peakCh';
cluster_channel = chanMap(cluster_channel);
clusterScore = clusterScore;

save([ses_name, '_Phy_', region], 'spikes_phy', 'clusterId', 'cluster_channel', 'clusterScore')
disp('Done')
end