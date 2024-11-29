function RawNeural_Data_mV = doFilters(RawNeural_Data_mV, fs)

disp('Applying Filters')

wo = 1990/(fs/2);
bw = wo/10;
[b1,a1] = iirnotch(wo,bw);
wo = 2000/(fs/2);
bw = wo/10;
[b2,a2] = iirnotch(wo,bw);
wo = 2010/(fs/2);
bw = wo/10;
[b3,a3] = iirnotch(wo,bw);
wo = 4000/(fs/2);
bw = wo/10;
[b4,a4] = iirnotch(wo,bw);
wo = 4020/(fs/2);
bw = wo/10;
[b5,a5] = iirnotch(wo,bw);
wo = 3080/(fs/2);
bw = wo/10;
[b6,a6] = iirnotch(wo,bw);

hpFilt = designfilt('highpassiir','FilterOrder',8, ...
    'PassbandFrequency',150,'PassbandRipple',0.2, ...
    'SampleRate',30000);

d1 = designfilt('bandstopiir', ...
    'PassbandFrequency1',1950,'StopbandFrequency1',1990,...
    'StopbandAttenuation',60,'PassbandRipple1',1, ...
    'PassbandRipple2', 1, 'StopbandFrequency2', 2010, ...
    'PassbandFrequency2', 2050, ...
    'DesignMethod','butter','SampleRate',fs);

d2 = designfilt('bandstopiir', ...
    'PassbandFrequency1',3950,'StopbandFrequency1',3990,...
    'StopbandAttenuation',60,'PassbandRipple1',1, ...
    'PassbandRipple2', 1, 'StopbandFrequency2', 4010, ...
    'PassbandFrequency2', 4050, ...
    'DesignMethod','butter','SampleRate',fs);

parpool; % Start the parallel pool
parfor i = 1:size(RawNeural_Data_mV, 1)
    display(['Filtering Ch# ',num2str(i)])
    % tic
    tmp = filtfilt(d1, RawNeural_Data_mV(i, :));
    tmp = filtfilt(d2, tmp);
    tmp = bandpass(tmp,[700 5000],fs)
    filtered_data6(i, :) = tmp;
end
delete(gcp('nocreate'))

RawNeural_Data_mV = filtered_data6;

if(saveFilt); display('Saving filtered data'); save(['Z_', sesName, '_', region, '_Filtered'], 'RawNeural_Data_mV','-v7.3'); end

end