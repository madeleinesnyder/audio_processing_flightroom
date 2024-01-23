%% Plot the X ms around a given "audio event" WITH THE FLIGHT DATA

% List of good segments (desired events)

% For a given bat, and session
Bat = 32626; logger = 15;
batdate = 221126;

save_audio_and_figure = 0;

% Audio and flight sampling rates
motu_Fs = 192000;
ciholas_Fs = 120;
mic_num = 4;
context_interval_ = 5000000;

%% Load in the audio data
audio_base = strcat('C:\Users\YartsevLabComputer5\Desktop\audio_data\',num2str(batdate),'\');
load(strcat(audio_base,'event_timestamps.mat'),'event_locs');

% Load in ttl data
alignment_='start';
load(strcat(audio_base,'ttl_locs.mat'));
load(strcat(audio_base,'ttl_first_sample.mat'));
load(strcat(audio_base,'ttl_last_sample.mat'));

%sorted_events = sorted_events;%-first_ttl_sample;
sorted_events = unique(sort([event_locs{1},event_locs{2},event_locs{3},event_locs{4}]));
%sorted_events = unique(sort([event_locs{1}]));

% Plot the sorted events and user-select an audio event timestamp to plot
% around.
figure(); hold on; plot(sorted_events)
title('Click on the plot to get which audio sample you want to examine');
% Wait for a mouse click and get the coordinates
[xClick, yClick] = ginput(1);
hold off;
desired_event = find(abs(sorted_events-yClick) == min(abs(sorted_events-yClick)));
desired_event = sorted_events(desired_event);

%desired_event = first_ttl_sample + 10;
%desired_event = 935703747;
%desired_event = round(desired_event+(8.41858*192000)+(1.614*192000));
%desired_event = round(desired_event+(5.59508*192000));
%desired_event = 1378936894;
%desired_event = 208464052;

% Sample
% Find which chunk and seg this event is from.
file_parts = dir(audio_base); 
file_parts = strsplit(file_parts(end-4).name,'_');
max_chunks = str2num(file_parts{end-2});
running_chunk_start = 0;
desired_chunk = []; chunk_start = []; amount_into_chunk = [];
for tt = 1:max_chunks

    chunk_filename = dir(strcat(audio_base,'/audioConCat_mic_1_segment_chunk_',num2str(tt),'_break_*.mat'));
    load(strcat(audio_base,chunk_filename.name));
    len_chunk = length(audioConCat);
    
    % Add the length of the chunk of audio to the total amount of length in the chunks so far
    running_chunk_start = running_chunk_start + len_chunk;

    if (desired_event < running_chunk_start) & (desired_event > (running_chunk_start-len_chunk))
        desired_chunk = tt;
        chunk_start = running_chunk_start-len_chunk;
        amount_into_chunk = desired_event - (running_chunk_start-len_chunk);
        break;
    end
end

% Plot an example
chunk = desired_chunk;  seg_start = amount_into_chunk; seg_end = amount_into_chunk+context_interval_;
event_vec_start = desired_event;
event_vec_end = desired_event+context_interval_;

%% Get filtered data for segment from all mics 
whitenedSignal_mm = {};
for mm = 1:mic_num 
    
    audio_file = dir(strcat(audio_base,'/audioConCat_mic_',num2str(mm),'_segment_chunk_',num2str(chunk),'_break_*.mat'));
    load(strcat(audio_base,audio_file.name));

    % Apply segment bounds to desired segment to plot 
    if seg_end>len_chunk
        signal = audioConCat(seg_start:len_chunk);
    else
        signal = audioConCat(seg_start:seg_end);
    end
    
    %% Remove 60 and 90 Hz noise from segment
    
    % Design a notch filter to remove 60 Hz noise
    wo_60 = 60/(motu_Fs/2);  % 60 Hz frequency in normalized frequency units
    bw_60 = wo_60/5;      % Bandwidth
    [b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
    filteredAudio_1 = filter(b, a, signal);
        
    % Notch filter for 90 Hz
    wo_90 = 90/(motu_Fs/2);  % 60 Hz frequency in normalized frequency units
    bw_90 = wo_90/5;      % Bandwidth
    [b, a] = iirnotch(wo_90, bw_90);  % IIR Notch filter design
    % Apply the notch filter to the audio data
    filteredAudio_2 = filter(b, a, filteredAudio_1);
    
    %% Whiten the signal segment (flattening the power spectrum)
    
    % FFT the signal and 
    fftSignal = fft(filteredAudio_2);
    powerSpectrum = abs(fftSignal).^2;
    
    whitenedFFTSignal = fftSignal ./ sqrt(powerSpectrum + eps); % eps is added for numerical stability
    whitenedSignal = real(ifft(whitenedFFTSignal));
    
    whitenedSignal_mm{mm} = whitenedSignal / max(abs(whitenedSignal));

end

whitenedSignal_allmics = cell2mat(whitenedSignal_mm);
sum_whitenedSignal_allmics = sum(whitenedSignal_allmics,2);
    
% above_zero_events = (sorted_events-event_vec_start);
% above_zero_events = above_zero_events(above_zero_events>0);
% above_zero_events = above_zero_events(above_zero_events<length(sum_whitenedSignal_allmics));
% above_zero_events_aligned_to_start_of_session = above_zero_events+event_vec_start;

%% Align the above_zero_events to the flight data 
above_zero_events = (sorted_events-event_vec_start); %above_zero_events = (scaled_zerod_loc_vec-event_vec_start);
above_zero_events = above_zero_events(above_zero_events>=0);
above_zero_events(above_zero_events==0) = [];
above_zero_events = above_zero_events(above_zero_events<length(sum_whitenedSignal_allmics));
%above_zero_events_aligned_to_start_of_session = above_zero_events+event_vec_start;
%[above_zero_events_aligned,ttl_events_aligned] = HumanBat_align_audio_to_ciholas(Bat,batdate,logger,alignment_,above_zero_events_aligned_to_start_of_session,pk_loc_vec,first_ttl_sample,last_ttl_sample);

% above_zero_events_aligned are the motu timestamps of the audio events in
% the chosen segment ALIGNED to the "alignement_" of the session.

%start_sample_ciholas = (event_vec_start-first_ttl_sample)/motu_Fs*ciholas_Fs;
%end_sample_ciholas = (event_vec_end-first_ttl_sample)/motu_Fs*ciholas_Fs;

start_sample_ciholas = round((event_vec_start-first_ttl_sample)/motu_Fs*ciholas_Fs);
if start_sample_ciholas == 0; start_sample_ciholas = 1; end;
end_sample_ciholas = round((event_vec_end-first_ttl_sample)/motu_Fs*ciholas_Fs);

%% Load in the flight data, extracted behavior, and ephys data
exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/');
    
% Load in the flight data 
load(strcat(exp_data_path,'ciholas/aligned_bat_position_data_',num2str(Bat),'.mat'));

% Load in the Extracted Behavior data
load(strcat(exp_data_path,'ciholas/Extracted_Behavior_', num2str(batdate),'_',num2str(Bat),'.mat'));

% Load in spike data
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));

%% Find the ciholas and ephys timestamps for this audio segment

figure(); 
subplot(1,2,1); hold on; 
set(gcf, 'Color', 'k');
xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);

time_vec = [1:length(ciholas_r(start_sample_ciholas:end_sample_ciholas,1))];
t_normalized = (time_vec - min(time_vec)) / (max(time_vec) - min(time_vec));

% Create a custom color map from light green to light orange
green = [144, 238, 144] / 255; % Light green in RGB
orange = [255, 165, 0] / 255; % Light orange in RGB
colorMap = [linspace(green(1), orange(1), 256)', ...
            linspace(green(2), orange(2), 256)', ...
            linspace(green(3), orange(3), 256)'];
colormap(colorMap);
xlabel("mm",'FontWeight','bold');
ylabel("mm",'FontWeight','bold');
title("32626 flight paths during audio segment",'FontWeight','bold','Color','w');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');

scatter3(ciholas_r(start_sample_ciholas:end_sample_ciholas,1),...
    ciholas_r(start_sample_ciholas:end_sample_ciholas,2),...
    ciholas_r(start_sample_ciholas:end_sample_ciholas,3), 10 ,t_normalized,'filled');
grid off;
hold off;

subplot(1,2,2);  hold on; 
set(gcf, 'Color', 'k');
xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);

time_vec = [1:length(ciholas_r(start_sample_ciholas:end_sample_ciholas,1))];
t_normalized = (time_vec - min(time_vec)) / (max(time_vec) - min(time_vec));

% Create a custom color map from light green to light orange
green = [144, 238, 144] / 255; % Light green in RGB
orange = [255, 165, 0] / 255; % Light orange in RGB
colorMap = [linspace(green(1), orange(1), 256)', ...
            linspace(green(2), orange(2), 256)', ...
            linspace(green(3), orange(3), 256)'];
colormap(colorMap);
title("14643 flight paths during audio segment",'FontWeight','bold','Color','w');
xlabel("mm",'FontWeight','bold');
ylabel("mm",'FontWeight','bold');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');

ciholas_r_obat = load(strcat(exp_data_path,'ciholas/aligned_bat_position_data_14643.mat'));

scatter3(ciholas_r_obat.ciholas_r(start_sample_ciholas:end_sample_ciholas,1),...
    ciholas_r_obat.ciholas_r(start_sample_ciholas:end_sample_ciholas,2),...
    ciholas_r_obat.ciholas_r(start_sample_ciholas:end_sample_ciholas,3), 10 ,t_normalized,'filled');
grid off;
hold off;

close all;

%% Plot the signal segment and the echolocations from each microphone

figure(); 
subplot(4,1,1); hold on; 
title("All mic echolocations plotted on filtered audio",'FontWeight','bold','Color','w');
for mm=1:mic_num
    plot(whitenedSignal_mm{mm});
end
ee = NaN(1,length(sum_whitenedSignal_allmics));
ee(above_zero_events) = 1;
%ee_spec = NaN(1,length(sum_whitenedSignal_allmics));
%ee_spec( desired_event-event_vec_start) = 1;
scatter([1:length(sum_whitenedSignal_allmics)],ee,10,'r');
%scatter([1:length(sum_whitenedSignal_allmics)],ee_spec,10,'filled','g');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gcf, 'Color', 'k'); 
xlabel("time (samples)",'FontWeight','bold');
ylabel("amplitude",'FontWeight','bold');
axis tight;
hold off;

% Below from fb_pretty_sonogram code in FinchScope repo
overlap=2000;
tscale=2;
N=2048;
nfft=2^nextpow2(N);
low=2.9;
high=10;

t_=-N/2+1:N/2;
sigma=(tscale/1e3)*motu_Fs;
w = exp(-(t_/sigma).^2);
dw = -2*w.*(t_/(sigma^2));
    
% Calculate spectrogram
[S,F,T,P]=spectrogram(sum_whitenedSignal_allmics,w,overlap,nfft,motu_Fs);
[S2]=spectrogram(sum_whitenedSignal_allmics,dw,overlap,nfft,motu_Fs);
IMAGE=100*((abs(S)+abs(S2))/2);
IMAGE=log(IMAGE+2);
IMAGE(IMAGE>high)=high;
IMAGE(IMAGE<low)=low;
IMAGE=IMAGE-low;
IMAGE=IMAGE/(high-low);
IMAGE=63*(IMAGE); % 63 is the cutoff for GIF
IMAGE_MODS = log(abs(IMAGE)+1e+2);
% Above from fb_pretty_sonogram code in FinchScope repo
IMAGE_MODS_SUM = sum(IMAGE_MODS);

spectrogram_ee = NaN(length(T),1);
spectrogram_ee_idxs = floor(above_zero_events/(length(sum_whitenedSignal_allmics)/length(T)));
if ~isempty(above_zero_events)
    if spectrogram_ee_idxs(1) == 0; spectrogram_ee_idxs(1) = 1; end;
end
spectrogram_ee(spectrogram_ee_idxs) = 60000;

% Plot the spectrogram with the peaks highlighted
subplot(4,1,2); hold on; 
colormap(hot)
imagesc(T,F,IMAGE_MODS); set(gca,'YDir','normal');
axis tight;
plot(T, spectrogram_ee','y*','MarkerSize',5); % Mark the peaks
ylabel('Frequency kHz','FontWeight','bold');
xlabel('time (s)','FontWeight','bold');
title('Spectrogram of whitened data','FontWeight','bold','Color','w');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gcf, 'Color', 'k'); 
hold off;

% Plot the wingbeat alignment
subplot(4,1,3); hold on; 

% Load in the Extracted Behavior data
load(strcat(exp_data_path,'ciholas/Extracted_Behavior_', num2str(batdate),'_',num2str(Bat),'.mat'));

% Load in spike data
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));

% First ephys TTL time
usec_of_first_ephys_TTL = (B_ephys_data.TT_unit(1).Timestamps(1)-B_ephys_data.TT_unit(1).AlignedTimestamps(1));
% ^^ This is the timestamp usec you want to subtract from all values in the CSC file to get the same alignment to the ephys spike file

% Load in the CSC.mat data
csc_files = dir(fullfile(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/*CSC*')));
load(strcat(csc_files(1).folder,'/',csc_files(1).name));
Fs = Estimated_channelFS_Transceiver(1);
raw_V = double(AD_count_int16*AD_count_to_uV_factor);
tstamps = Timestamps_of_first_samples_usec(1)+[0:length(raw_V)-1]/Fs*1e6;
tstamps_shift = tstamps-usec_of_first_ephys_TTL;   
t_stamps_g0 = find(tstamps_shift > 0); t_g0 = t_stamps_g0(1); t_stamps_g0=[];
tstamps = []; 
raw_V_shift = raw_V(t_g0:end); raw_V = [];

% % Plot the ciholas acceleration value for seg_amt samples ciholas seconds is /120
% % ephys seconds is /1e6
% seg_amt_sec = 780; % Take first 780 seconds of the session
% ciholas_sample_amt = seg_amt_sec*120;
% ephys_sample_amt = seg_amt_sec*1e6;

% Downsample the ephys data
t_sample_length = t(end)/length(t);
extra_time = tstamps_shift(end)/1e6-t(end);
extra_samples = round(extra_time/t_sample_length);
desired_ephys_trail_sample_length = length(t)+extra_samples;

oL = length(raw_V_shift);
y = interp1(1:oL,raw_V_shift,linspace(1,oL,desired_ephys_trail_sample_length)); y = y/3000;

% Get indexes to start and stop the plots
t_start_sec = start_sample_ciholas/ciholas_Fs;
t_end_sec = end_sample_ciholas/ciholas_Fs;
t_start_index = find(abs(t-t_start_sec) == min(abs(t-t_start_sec)));
t_end_index = find(abs(t-t_end_sec) == min(abs(t-t_end_sec)));

% Plot ciholas data, wingbeats, and ephys!
[up,lo] = envelope(a_flt+normrnd(0,1e-3,length(a_flt),1),ciholas_Fs/10,'peak');     % Envelope of the acceleration signal (noise addedd to avoid problems with splining)
env = normalize(up - lo,'range');                                           % Amplitude of the envelope
env_th = otsuthresh(histcounts(env));                                       % Threshold (based on Otsu method). Can be set at 0.35
wBeats = movsum(env>env_th,1.67*ciholas_Fs)>ciholas_Fs/6; 
area(t(t_start_index:t_end_index),wBeats(t_start_index:t_end_index)*3,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(t(t_start_index:t_end_index),normalize(v_abs(t_start_index:t_end_index),'range'));
plot(t(t_start_index:t_end_index),r(t_start_index:t_end_index,1),t(t_start_index:t_end_index),r(t_start_index:t_end_index,2));
plot(t(t_start_index:t_end_index),normalize(a_flt(t_start_index:t_end_index),'range',[-1 1]));
plot(t(t_start_index:t_end_index),normalize(movsum(env(t_start_index:t_end_index)>env_th,1.67*ciholas_Fs),'range'));
plot(t(t_start_index:t_end_index),up(t_start_index:t_end_index));
plot(t(t_start_index:t_end_index),lo(t_start_index:t_end_index));
plot(t(t_start_index:t_end_index),y(1:length(t(t_start_index:t_end_index))));
Z = round(above_zero_events/motu_Fs*ciholas_Fs); Z(Z==0)=1;
binZ = NaN(1,length([t_start_index:t_end_index])); 
binZ(Z) = 4;
scatter(t(t_start_index:t_end_index),binZ,'*r');
[y_up,y_lo] = envelope(y(1:length(t(t_start_index:t_end_index))),ciholas_Fs/10,'peak');     % Envelope of the acceleration signal (noise addedd to avoid problems with splining)
y_env = normalize(y_up - y_lo,'range'); 
plot(t(t_start_index:t_end_index),y_up);
plot(t(t_start_index:t_end_index),y_lo);
xlabel("time (samples)",'FontWeight','bold');
ylabel("amplitude",'FontWeight','bold');
title("Bat 32626 wingbeat sync",'FontWeight','bold','Color','w');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gcf, 'Color', 'k'); 
axis tight;
hold off;

clear a_flt wBeats up lo env 

% Plot the wingbeat alignment
subplot(4,1,4); hold on; 

obat = 14643; 
ologger = 13;

% Load in the Extracted Behavior data
load(strcat(exp_data_path,'ciholas/Extracted_Behavior_', num2str(batdate),'_',num2str(obat),'.mat'));

% Load in spike data
load(strcat(exp_data_path,'ephys/logger',num2str(ologger),'/extracted_data/B_ephys_data_aligned.mat'));

% First ephys TTL time
usec_of_first_ephys_TTL = (B_ephys_data.TT_unit(1).Timestamps(1)-B_ephys_data.TT_unit(1).AlignedTimestamps(1));
% ^^ This is the timestamp usec you want to subtract from all values in the CSC file to get the same alignment to the ephys spike file

% Load in the CSC.mat data
csc_files = dir(fullfile(strcat(exp_data_path,'ephys/logger',num2str(ologger),'/extracted_data/*CSC*')));
load(strcat(csc_files(1).folder,'/',csc_files(1).name));
Fs = Estimated_channelFS_Transceiver(1);
raw_V = double(AD_count_int16*AD_count_to_uV_factor);
tstamps = Timestamps_of_first_samples_usec(1)+[0:length(raw_V)-1]/Fs*1e6;
tstamps_shift = tstamps-usec_of_first_ephys_TTL;   
t_stamps_g0 = find(tstamps_shift > 0); t_g0 = t_stamps_g0(1); t_stamps_g0=[];
tstamps = []; 
raw_V_shift = raw_V(t_g0:end); raw_V = [];

% Downsample the ephys data
t_sample_length = t(end)/length(t);
extra_time = tstamps_shift(end)/1e6-t(end);
extra_samples = round(extra_time/t_sample_length);
desired_ephys_trail_sample_length = length(t)+extra_samples;

oL = length(raw_V_shift);
y = interp1(1:oL,raw_V_shift,linspace(1,oL,desired_ephys_trail_sample_length)); y = y/3000;

% Get indexes to start and stop the plots
t_start_sec = (start_sample_ciholas)/ciholas_Fs;
t_end_sec = (end_sample_ciholas)/ciholas_Fs;
t_start_index = find(abs(t-t_start_sec) == min(abs(t-t_start_sec)));
t_end_index = find(abs(t-t_end_sec) == min(abs(t-t_end_sec)));

% Plot ciholas data, wingbeats, and ephys!
[up,lo] = envelope(a_flt+normrnd(0,1e-3,length(a_flt),1),ciholas_Fs/10,'peak');     % Envelope of the acceleration signal (noise addedd to avoid problems with splining)
env = normalize(up - lo,'range');                                           % Amplitude of the envelope
env_th = otsuthresh(histcounts(env));                                       % Threshold (based on Otsu method). Can be set at 0.35
wBeats = movsum(env>env_th,1.67*ciholas_Fs)>ciholas_Fs/6; 
area(t(t_start_index:t_end_index),wBeats(t_start_index:t_end_index)*3,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(t(t_start_index:t_end_index),normalize(v_abs(t_start_index:t_end_index),'range'));
plot(t(t_start_index:t_end_index),r(t_start_index:t_end_index,1),t(t_start_index:t_end_index),r(t_start_index:t_end_index,2));
plot(t(t_start_index:t_end_index),normalize(a_flt(t_start_index:t_end_index),'range',[-1 1]));
plot(t(t_start_index:t_end_index),normalize(movsum(env(t_start_index:t_end_index)>env_th,1.67*ciholas_Fs),'range'));
plot(t(t_start_index:t_end_index),up(t_start_index:t_end_index));
plot(t(t_start_index:t_end_index),lo(t_start_index:t_end_index));
plot(t(t_start_index:t_end_index),y(1:length(t(t_start_index:t_end_index))));
Z = round(above_zero_events/motu_Fs*ciholas_Fs); Z(Z==0)=1;
binZ = NaN(1,length([t_start_index:t_end_index])); 
binZ(Z) = 4;
scatter(t(t_start_index:t_end_index),binZ,'*r');
[y_up,y_lo] = envelope(y(1:length(t(t_start_index:t_end_index))),ciholas_Fs/10,'peak');     % Envelope of the acceleration signal (noise addedd to avoid problems with splining)
y_env = normalize(y_up - y_lo,'range'); 
plot(t(t_start_index:t_end_index),y_up);
plot(t(t_start_index:t_end_index),y_lo);
xlabel("time (samples)",'FontWeight','bold');
ylabel("amplitude",'FontWeight','bold');
title("Bat 14643 wingbeat sync",'FontWeight','bold','Color','w');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gcf, 'Color', 'k'); 
axis tight;
hold off;

% Clip that segment out and put it on the server 

% if save_audio_and_figure == 1
% 
%     clip_name = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/audio_clips/',num2str(batdate),'_',num2str(Bat),'_event_start_sample_',num2str(desired_event),'.wav');
%     audiowrite(clip_name, sum_whitenedSignal_allmics, motu_Fs);
%     
%     clip_name_figure = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/audio_clips/',num2str(batdate),'_',num2str(Bat),'_event_start_sample_',num2str(desired_event),'.svg');
%     saveas(gcf, clip_name_figure, 'svg');
% end








