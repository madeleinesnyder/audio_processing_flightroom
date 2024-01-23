%% Load in a chunk of audio and save the timestamps of echolocations to a file

batdate = 221125;
mic = 2;
chunk = 5;

audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');
load(strcat(audio_base,'/audioConCat_mic_',num2str(mic),'_segment_chunk_',num2str(chunk),'_break_0.mat'));
Fs = 192000;

%% Remove 60 and 90 Hz noise
% Design a notch filter to remove 60 Hz noise
wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_60 = wo_60/5;      % Bandwidth
[b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
filteredAudio_1 = filter(b, a, audioConCat);

% Notch filter for 90 Hz
wo_90 = 90/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_90 = wo_90/5;      % Bandwidth
[b, a] = iirnotch(wo_90, bw_90);  % IIR Notch filter design
% Apply the notch filter to the audio data
filteredAudio_2 = filter(b, a, filteredAudio_1);

%% Whiten the signal (flattening the power spectrum)

% FFT the signal and 
fftSignal = fft(filteredAudio_2);
powerSpectrum = abs(fftSignal).^2;

whitenedFFTSignal = fftSignal ./ sqrt(powerSpectrum + eps); % eps is added for numerical stability
whitenedSignal = real(ifft(whitenedFFTSignal));

whitenedSignal = whitenedSignal / max(abs(whitenedSignal));

%sound(whitenedSignal, Fs); % To listen

%% Bandpass filter the segment to include only from 4000 to 8000 Hz

% Bandpass between 2Hz and 6kHz
% Design a bandpass filter
d = designfilt('bandpassiir', ...
               'FilterOrder', 6, ...         % Filter order
               'HalfPowerFrequency1', 4000, ... % Lower cutoff frequency
               'HalfPowerFrequency2', 8000, ... % Upper cutoff frequency
               'SampleRate', Fs);              % Sampling rate

bp_echo_whitenedSignal = filter(d, whitenedSignal);

%% FFT visualization

addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/shared_utils/FinchScope'));

echo_segment_start = 135000000;
echo_segment_end = 140000000;
% Sample signal
t = 0:1/Fs:1; % Time vector (1 second)
signal = bp_echo_whitenedSignal(echo_segment_start:echo_segment_end); % Sine wave
signal_no_echo_bp = whitenedSignal(echo_segment_start:echo_segment_end); % Sine wave
save('sample_chunk_5_signal.mat','signal');
save('sample_chunk_5_signal_no_echo_bp.mat','signal_no_echo_bp');

% % Generate Spectrogram of the bandpass filtered and whitened segment
% figure('name', 'Raw Whitened and 4-8kHz bp data'); 
% [IMAGE,F,T] = fb_pretty_sonogram(signal./abs(max(signal)),Fs,'low',2.9,'zeropad',0);
% colormap(hot)
% imagesc(T,F,log(abs(IMAGE)+1e+2));set(gca,'YDir','normal');
% ylabel('kHz')
% xlabel('time (s)');
% title('Spectrogram of manually selected audio subsegment');
% 
% % Generate Spectrogram of the whitened segment 
% figure('name', 'Raw Whitened Data'); 
% [IMAGE,F,T] = fb_pretty_sonogram(signal_no_echo_bp./abs(max(signal_no_echo_bp)),Fs,'low',2.9,'zeropad',0);
% colormap(hot)
% imagesc(T,F,log(abs(IMAGE)+1e+2));set(gca,'YDir','normal');
% ylabel('kHz')
% xlabel('time (s)');
% title('Spectrogram of manually selected audio subsegment');

%% Extract peaks from the bandpass and whitened signal 

% There is a time-frequency trade-off. The longer the window, the higher the
% frequency resoltuion but the poorer the time resolution. 

% Below from fb_pretty_sonogram code in FinchScope repo
overlap=2000;
tscale=2;
N=2048;
nfft=2^nextpow2(N);
low=2.9;
high=10;

t=-N/2+1:N/2;
sigma=(tscale/1e3)*Fs;
w = exp(-(t/sigma).^2);
dw = -2*w.*(t/(sigma^2));

% Calculate spectrogram
[S,F,T,P]=spectrogram(signal_no_echo_bp,w,overlap,nfft,Fs);
[S2]=spectrogram(signal_no_echo_bp,dw,overlap,nfft,Fs);
IMAGE=100*((abs(S)+abs(S2))/2);
IMAGE=log(IMAGE+2);
IMAGE(IMAGE>high)=high;
IMAGE(IMAGE<low)=low;
IMAGE=IMAGE-low;
IMAGE=IMAGE/(high-low);
IMAGE=63*(IMAGE); % 63 is the cutoff for GIF
IMAGE_MODS = log(abs(IMAGE)+1e+2);
% Above from fb_pretty_sonogram code in FinchScope repo

% S: Short-time Fourier Transform
% F: Frequency vector
% T: Time vector
% P: Power spectral density

%% Get the peaks of the summed frequency bands

IMAGE_MODS_SUM = sum(IMAGE_MODS);
threshold = 1;
peak_dist_threshold_seconds = 0.0012*Fs;
peakFreqs = NaN(size(P, 1), 1);
[idx, loc] = findpeaks(abs(IMAGE_MODS_SUM), 'MinPeakProminence', threshold,'MinPeakDistance',peak_dist_threshold_seconds);

% Turn locations into a binary vector and scale the 1's to high frequency
% for plotting purposes.
peakFreqs(find(isnan(loc))) = [];
pf_binary = zeros(length(T),1); pf_binary(loc) = 1;
pf_binary(find(pf_binary==1)) = 63000;

% Save the times of the echolocation peaks in the long vector. 
scaled_loc = round(loc*(length(signal_no_echo_bp)/length(T)));     % Scale the locations of the echolocations according to the samples in the real data vector

scaled_loc_binary = NaN(length(signal_no_echo_bp),1); scaled_loc_binary(scaled_loc) = 0.3;
figure(); hold on; plot(signal_no_echo_bp); 
scatter(1:length(signal_no_echo_bp),scaled_loc_binary);

% Plot the spectrogram with the peaks highlighted
figure(); 
subplot(2,1,1); hold on; 
colormap(hot)
imagesc(T,F,IMAGE_MODS); set(gca,'YDir','normal');
plot(T, pf_binary, 'y*', 'MarkerSize', 1); % Mark the peaks
ylabel('Frequency kHz')
xlabel('time (s)');
title('Spectrogram of whitened data');
hold off;

subplot(2,1,2); hold on;
title("Whitened signal with echolocation peaks");
xlabel("Time (samples) at 192000Hz");
ylabel("Amplitude");
plot(signal_no_echo_bp); 
scatter(1:length(signal_no_echo_bp),scaled_loc_binary);






