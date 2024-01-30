%% Script to bandpass, downsample, and then event detect

%% For a given bat, and session
Bat = 32626; logger = 15; batdate = 221126;

desired_event = 208464052;

save_audio_and_figure = 0;

% Audio and flight sampling rates
motu_Fs = 192000;
ciholas_Fs = 120;
mic_num = 4;
context_interval_ = 5000000;

% Spectrogram hyperparameters 
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

%% Load in the audio data
audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');
load(strcat(audio_base,'event_timestamps.mat'),'event_locs');
load(strcat(audio_base,'ttl_first_sample.mat'));

% Load in ttl data
alignment_='start';
load(strcat(audio_base,'ttl_locs.mat'));
load(strcat(audio_base,'ttl_first_sample.mat'));
load(strcat(audio_base,'ttl_last_sample.mat'));

% Find which chunk and seg this event is from.
file_parts = dir(audio_base); 
file_parts = strsplit(file_parts(end-3).name,'_');
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
    
    %% Echolocation detection processing: 
    % bandpass from 20kHz to 50kHz
    fpass_lo = 40000;
    fpass_hi = 60000;
    db_noise = 80;
    echo_amp_env_cutoff = 0.001;

    echoFilter = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',fpass_lo,'CutoffFrequency2',fpass_hi,'SampleRate',motu_Fs);
    signal_echobp = filtfilt(echoFilter,signal);

    % Make the spectrogram of the data (Input the bandpass filtered signal)
    [to, fo, logB, pg, tError, fError] = HumanBat_Julie_Spec_Only_Bats(signal_echobp, motu_Fs, db_noise,fpass_hi);

    % Run RMS amplitude calcuation (Input the bandpass filtered signal)
    FigFlag = 1;
    FS_high_cutoff = 50; % Play with this, set higher or lower depending on kind of input. Maybe 50? 100? for echolocations
    FS_env = 100; % Default
    [Amp_env_voltage, Power_env] = HumanBat_running_rms(signal_echobp, motu_Fs, FS_high_cutoff, FS_env, 'filter', FigFlag);
    Amp_env_voltage = Amp_env_voltage*10;
    t2 = (0:(length(Amp_env_voltage)-1))*motu_Fs/(FS_env);

    % 2. Use findpeaks to find times of echocoloations 
    mpp = 0.0008;
    mph = echo_amp_env_cutoff;
    mpd = ceil((192000*0.02)/(length(signal_echobp)/length(t2)));
    [pks,locs] = findpeaks(whitenedSignal_final,'MinPeakProminence',mpp,'MinPeakHeight',mph,'MinPeakDistance',mpd)
    Amp_env_binary2 = NaN(1,length(signal_echobp)); echolocation_timestamps = round(locs*(length(signal_echobp)/length(whitenedSignal_final))); Amp_env_binary2(echolocation_timestamps) = 0;
    
    % Plot
    figure(); hold on;
    plot(whitenedSignal_final);
    scatter([1:5000001],Amp_env_binary2,5,'filled','red');

    % echolocation_timestamps contains the timestamps of detected echolocations!

    %% Feeder detection processing: 
    % bandpass from 2500Hz to 20kHz
    fpass_lo = 2500;
    fpass_hi = 10000;
    db_noise = 80;
    feeder_amp_env_cutoff = 0.3;

    feederFilter = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',fpass_lo,'CutoffFrequency2',fpass_hi,'SampleRate',motu_Fs);
    signal_feederbp = filtfilt(feederFilter,signal);

    % Make the spectrogram of the data (Input the bandpass filtered signal)
    [to, fo, logB, pg, tError, fError] = HumanBat_Julie_Spec_Only_Bats(signal_feederbp, motu_Fs, db_noise, fpass_hi);

    % Run RMS amplitude calcuation (Input the bandpass filtered signal)
    FigFlag = 1;
    FS_high_cutoff = 100; % Play with this, set higher or lower depending on kind of input. Maybe 50? 100? for echolocations
    FS_env = 100; % Default
    [Amp_env_voltage, Power_env] = HumanBat_running_rms(signal_feederbp, motu_Fs, FS_high_cutoff, FS_env, 'filter', FigFlag);
    Amp_env_voltage = Amp_env_voltage;
    
    mpp = 0.2;
    mph = feeder_amp_env_cutoff;
    mpd = round(192000/200);
    [pks,locs] = findpeaks(Amp_env_voltage,'MinPeakProminence',mpp,'MinPeakHeight',mph,'MinPeakDistance',mpd);
    Amp_env_binary2 = NaN(1,length(signal_feederbp)); feeder_timestamps = round(locs*(length(signal_feederbp)/length(Amp_env_voltage))); Amp_env_binary2(feeder_timestamps) = 0;
    
    % Plot
    figure(); hold on;
    t2 = (0:(length(Amp_env_voltage)-1))*motu_Fs/(FS_env);
    plot(signal_feederbp);
    scatter([1:5000001],Amp_env_binary2,5,'filled','red');

    % feeder_timestamps contains timestamps of feeder clicks!
   
    %% Human click detection processing: 
    % bandpass from 10kHz to 40kHz
    fpass_lo = 10000;
    fpass_hi = 40000;
    db_noise = 80;
    echo_amp_env_cutoff = 0.005;

    humanclickFilter = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',fpass_lo,'CutoffFrequency2',fpass_hi,'SampleRate',motu_Fs);
    signal_humanclickbp = filtfilt(humanclickFilter,signal);

    % Make the spectrogram of the data (Input the bandpass filtered signal)
    [to, fo, logB, pg, tError, fError] = HumanBat_Julie_Spec_Only_Bats(signal_humanclickbp, motu_Fs, db_noise,fpass_hi);

    % Run RMS amplitude calcuation (Input the bandpass filtered signal)
    FigFlag = 1;
    FS_high_cutoff = 50; % Play with this, set higher or lower depending on kind of input. Maybe 50? 100? for echolocations
    FS_env = 100; % Default
    [Amp_env_voltage, Power_env] = HumanBat_running_rms(signal_humanclickbp, motu_Fs, FS_high_cutoff, FS_env, 'filter', FigFlag);
    Amp_env_voltage = Amp_env_voltage*10;
    
    % 2. Use findpeaks to find times of echocoloations 
    mpp = 0.001;
    mph = echo_amp_env_cutoff;
    mpd = ceil((192000*0.02)/(length(signal_humanclickbp)/length(t2)));
    [pks,locs] = findpeaks(Amp_env_voltage,'MinPeakProminence',mpp,'MinPeakHeight',mph,'MinPeakDistance',mpd)
    Amp_env_binary2 = NaN(1,length(signal_humanclickbp)); humanclick_timestamps = round(locs*(length(signal_humanclickbp)/length(Amp_env_voltage))); Amp_env_binary2(humanclick_timestamps) = 0;
    
    % Plot
    figure(); hold on;
    t2 = (0:(length(Amp_env_voltage)-1))*motu_Fs/(FS_env);
    plot(signal_humanclickbp);
    scatter([1:5000001],Amp_env_binary2,5,'filled','red');

    % humanclick_timestamps contains the timestamps of detected human clicks!

    %% Human voice detection processing: 
    % Remove 60Hz noise from segment
    
    % Apply the Matlab function for speech detection.
    % Go to python script HumanBat_HumanSpeechRecognition.py

    %% Filter the signal segment for plotting.
    % Design a notch filter to remove 60 Hz noise
    wo_60 = 60/(motu_Fs/2);  % 60 Hz frequency in normalized frequency units
    bw_60 = wo_60/5;      % Bandwidth
    [b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
    filteredAudio_1 = filter(b, a, signal);
        
    % FFT the signal and whiten
    fftSignal = fft(filteredAudio_1);
    powerSpectrum = abs(fftSignal).^2;
    
    whitenedFFTSignal = fftSignal ./ sqrt(powerSpectrum + eps); % eps is added for numerical stability
    whitenedSignal = real(ifft(whitenedFFTSignal));
    
    whitenedSignal_mm{mm} = whitenedSignal / max(abs(whitenedSignal));
end

whitenedSignal_allmics = cell2mat(whitenedSignal_mm);
sum_whitenedSignal_allmics = sum(whitenedSignal_allmics,2);
 
%% Plot the events in different colors on the raw data and spectogram!

% Zero the timestamps detected above to the start of the signal vector
echo_events = echolocation_timestamps-event_vec_start; echo_events = echo_events(echo_events>=0); echo_events(echo_events==0) = []; echo_events = echo_events(echo_events<length(signal));
feeder_events = feeder_timestamps-event_vec_start; feeder_events = feeder_events(feeder_events>=0); feeder_events(feeder_events==0) = []; feeder_events = feeder_events(feeder_events<length(signal));
humanclick_events = humanclick_timestamps-event_vec_start; humanclick_events = humanclick_events(humanclick_events>=0); humanclick_events(humanclick_events==0) = []; humanclick_events = humanclick_events(humanclick_events<length(signal));

% Figure out where in ciholas we are
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

%% Plot the signal segment and the echolocations from each microphone

figure(); 
subplot(4,1,1); hold on; 
title("All mic echolocations plotted on filtered audio",'FontWeight','bold','Color','w');
for mm=1:mic_num
    plot(whitenedSignal_mm{mm});
end
ee_echo = NaN(1,length(sum_whitenedSignal_allmics)); ee_feeder = NaN(1,length(sum_whitenedSignal_allmics)); ee_humanclick = NaN(1,length(sum_whitenedSignal_allmics));
ee_echo(echo_events) = 1; ee_feeder(feeder_events) = 1; ee_humanclick(humanclick_events) = 1;
scatter([1:length(sum_whitenedSignal_allmics)],ee_echo,10,'g');
scatter([1:length(sum_whitenedSignal_allmics)],ee_humanclick,10,'y');
scatter([1:length(sum_whitenedSignal_allmics)],ee_feeder,10,'orange');
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

spectrogram_ee_echo = NaN(length(T),1); spectrogram_ee_feeder = NaN(length(T),1); spectrogram_ee_humanclick = NaN(length(T),1);
spectrogram_ee_echo_idxs = floor(echo_events/(length(sum_whitenedSignal_allmics)/length(T))); spectrogram_ee_feeder_idxs = floor(feeder_events/(length(sum_whitenedSignal_allmics)/length(T))); spectrogram_ee_hc_idxs = floor(humanclick_events/(length(sum_whitenedSignal_allmics)/length(T)));
if ~isempty(echo_events)
    if spectrogram_ee_echo_idxs(1) == 0; spectrogram_ee_echo_idxs(1) = 1; end;
end
if ~isempty(feeder_events)
    if spectrogram_ee_feeder_idxs(1) == 0; spectrogram_ee_feeder_idxs(1) = 1; end;
end
if ~isempty(humanclick_events)
    if spectrogram_ee_hc_idxs(1) == 0; spectrogram_ee_hc_idxs(1) = 1; end;
end
spectrogram_ee_echo(spectrogram_ee_echo_idxs) = 60000; spectrogram_ee_feeder(spectrogram_ee_feeder_idxs) = 60000; spectrogram_ee_hc(spectrogram_ee_hc_idxs) = 60000;

% Plot the spectrogram with the peaks highlighted
subplot(4,1,2); hold on; 
colormap(hot)
imagesc(T,F,IMAGE_MODS); set(gca,'YDir','normal');
axis tight;
plot(T, spectrogram_ee_echo','g*','MarkerSize',5); % Mark the peaks
plot(T, spectrogram_ee_feeder','o*','MarkerSize',5); % Mark the peaks
plot(T, spectrogram_ee_hc','y*','MarkerSize',5); % Mark the peaks
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








