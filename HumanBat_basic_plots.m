%% Make basic plots for a segment of unfiltered audio data 

%% Load the data
batdate = 221126;
mic = 1;
chunk = 3;
chunk = 5;
Fs = 192000;
    
%% Load in the audio data
audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');
chunk_filename = dir(strcat(audio_base,'/audioConCat_mic_',num2str(mic),'_segment_chunk_',num2str(chunk),'_break_*.mat'));
load(strcat(audio_base,chunk_filename.name));

seg_start_in_chunk = 30000000;
seg_end_in_chunk = 40000000;

signal = audioConCat(seg_start_in_chunk:seg_end_in_chunk);

%% Apply filters and repeat the process

% Design a notch filter to remove 60 Hz noise
wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_60 = wo_60/5;      % Bandwidth
[b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
filteredAudio_1 = filter(b, a, signal);

% Notch filter for 90 Hz
wo_90 = 90/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_90 = wo_90/5;      % Bandwidth
[b, a] = iirnotch(wo_90, bw_90);  % IIR Notch filter design
filteredAudio_2 = filter(b, a, filteredAudio_1);

%% Whiten the signal whole chunk (flattening the power spectrum)

% FFT the signal and 
fftSignal = fft(filteredAudio_2);
powerSpectrum = abs(fftSignal).^2;
whitenedFFTSignal = fftSignal ./ sqrt(powerSpectrum + eps); % eps is added for numerical stability
whitenedSignal = real(ifft(whitenedFFTSignal));
whitenedSignal = whitenedSignal / max(abs(whitenedSignal));

% Calculate spectrogram of filtered data
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
[S,F,T,P]=spectrogram(whitenedSignal,w,overlap,nfft,Fs);
[S2]=spectrogram(whitenedSignal,dw,overlap,nfft,Fs);
IMAGE=100*((abs(S)+abs(S2))/2);
IMAGE=log(IMAGE+2);
IMAGE(IMAGE>high)=high;
IMAGE(IMAGE<low)=low;
IMAGE=IMAGE-low;
IMAGE=IMAGE/(high-low);
IMAGE=63*(IMAGE); % 63 is the cutoff for GIF
IMAGE_MODS = log(abs(IMAGE)+1e+2);
% Above from fb_pretty_sonogram code in FinchScope repo
clear IMAGE S P;

%% Plot the audio filtered
figure(); 
subplot(2,1,1); hold on; 
set(gcf, 'Color', 'k'); hold on;
colormap(hot);
plot(whitenedSignal);
set(gca,'YDir','normal');
axis tight;
ylabel('time (samples)','FontWeight','bold');
xlabel('amplitude','FontWeight','bold');
title('spectrogram of whitened data','FontWeight','bold');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');
hold off;

%% Plot spectrogram
        
subplot(2,1,2); hold on; 
set(gcf, 'Color', 'k'); hold on;
colormap(hot)
imagesc(T,F,IMAGE_MODS); set(gca,'YDir','normal');
axis tight;
ylabel('frequency kHz','FontWeight','bold');
xlabel('time (s)','FontWeight','bold');
title('spectrogram of whitened data','FontWeight','bold');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');
hold off;


%% Plot filtered and unfiltered for a specific segment

ssic = seg_start_in_chunk+18.45*Fs; % feeder click (chunk 3; batdate = 221126; mic = 1; seg_start_in_chunk = 30000000; seg_end_in_chunk = 40000000;
seic = seg_start_in_chunk+19.62*Fs; % feeder click ^^

%ssic = seg_start_in_chunk+45.5*Fs; % echolocation train (chunk 3; batdate = 221126; mic = 1; seg_start_in_chunk = 30000000; seg_end_in_chunk = 40000000;
%seic = seg_start_in_chunk+47.5*Fs; % echolocation train ^^

%ssic = 159797612-(192000*3.6); % human click train (chunk 5; batdate = 221126; mic = 1)
%seic = 159797612+192000*1; % human click train ^^

ssic = 160871867; % madeleine "Ceasar" (chunk 5; batdate = 221126; mic = 1)
seic = 160871867+192000*1.5; % madeleine "Ceasar"

signal_fragment = audioConCat(ssic:seic);

% Calculate spectrogram of filtered data
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
[S,F,T,P]=spectrogram(signal_fragment,w,overlap,nfft,Fs);
[S2]=spectrogram(signal_fragment,dw,overlap,nfft,Fs);
IMAGE=100*((abs(S)+abs(S2))/2);
IMAGE=log(IMAGE+2);
IMAGE(IMAGE>high)=high;
IMAGE(IMAGE<low)=low;
IMAGE=IMAGE-low;
IMAGE=IMAGE/(high-low);
IMAGE=63*(IMAGE); % 63 is the cutoff for GIF
IMAGE_MODS_NOFILTER = log(abs(IMAGE)+1e+2);
% Above from fb_pretty_sonogram code in FinchScope repo
clear IMAGE S P;

% Design a notch filter to remove 60 Hz noise
wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_60 = wo_60/5;      % Bandwidth
[b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
filteredAudio_1 = filter(b, a, signal_fragment);

% Notch filter for 90 Hz
wo_90 = 90/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_90 = wo_90/5;      % Bandwidth
[b, a] = iirnotch(wo_90, bw_90);  % IIR Notch filter design
filteredAudio_2 = filter(b, a, filteredAudio_1);

% Calculate spectrogram
[S,F,T,P]=spectrogram(filteredAudio_2,w,overlap,nfft,Fs);
[S2]=spectrogram(filteredAudio_2,dw,overlap,nfft,Fs);
IMAGE=100*((abs(S)+abs(S2))/2);
IMAGE=log(IMAGE+2);
IMAGE(IMAGE>high)=high;
IMAGE(IMAGE<low)=low;
IMAGE=IMAGE-low;
IMAGE=IMAGE/(high-low);
IMAGE=63*(IMAGE); % 63 is the cutoff for GIF
IMAGE_MODS_BP = log(abs(IMAGE)+1e+2);
% Above from fb_pretty_sonogram code in FinchScope repo
clear IMAGE S P;

%% Whiten the signal whole chunk (flattening the power spectrum)

% FFT the signal and 
fftSignal = fft(filteredAudio_2);
powerSpectrum = abs(fftSignal).^2;
whitenedFFTSignal = fftSignal ./ sqrt(powerSpectrum + eps); % eps is added for numerical stability
whitenedSignal = real(ifft(whitenedFFTSignal));
whitenedSignal = whitenedSignal / max(abs(whitenedSignal));

% Calculate spectrogram
[S,F,T,P]=spectrogram(whitenedSignal,w,overlap,nfft,Fs);
[S2]=spectrogram(whitenedSignal,dw,overlap,nfft,Fs);
IMAGE=100*((abs(S)+abs(S2))/2);
IMAGE=log(IMAGE+2);
IMAGE(IMAGE>high)=high;
IMAGE(IMAGE<low)=low;
IMAGE=IMAGE-low;
IMAGE=IMAGE/(high-low);
IMAGE=63*(IMAGE); % 63 is the cutoff for GIF
IMAGE_MODS_WHITE = log(abs(IMAGE)+1e+2);
% Above from fb_pretty_sonogram code in FinchScope repo
clear IMAGE S P;

%% Plot all on one 
figure(); 

set(gcf, 'Color', 'k'); 
subplot(2,3,1); hold on;
colormap(hot);
plot(signal_fragment);
set(gca,'YDir','normal');
axis tight;
xlabel('time (samples)','FontWeight','bold');
ylabel('amplitude','FontWeight','bold');
title('spectrogram of whitened data','FontWeight','bold');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');
hold off;

subplot(2,3,2); hold on;
colormap(hot);
plot(filteredAudio_2);
set(gca,'YDir','normal');
axis tight;
xlabel('time (samples)','FontWeight','bold');
ylabel('amplitude','FontWeight','bold');
title('spectrogram of whitened data','FontWeight','bold');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');
hold off;

subplot(2,3,3); hold on;
colormap(hot);
plot(whitenedSignal);
set(gca,'YDir','normal');
axis tight;
xlabel('time (samples)','FontWeight','bold');
ylabel('amplitude normalized','FontWeight','bold');
title('spectrogram of whitened data','FontWeight','bold');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');
hold off;

subplot(2,3,4); hold on;
colormap(hot)
imagesc(T,F,IMAGE_MODS_NOFILTER); set(gca,'YDir','normal');
axis tight;
ylabel('frequency kHz','FontWeight','bold');
xlabel('time (s)','FontWeight','bold');
title('spectrogram of whitened data','FontWeight','bold');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');
hold off;

subplot(2,3,5); hold on;
colormap(hot)
imagesc(T,F,IMAGE_MODS_BP); set(gca,'YDir','normal');
axis tight;
ylabel('frequency kHz','FontWeight','bold');
xlabel('time (s)','FontWeight','bold');
title('spectrogram of whitened data','FontWeight','bold');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');
hold off;

subplot(2,3,6); hold on;
colormap(hot)
imagesc(T,F,IMAGE_MODS_WHITE); set(gca,'YDir','normal');
axis tight;
ylabel('frequency kHz','FontWeight','bold');
xlabel('time (s)','FontWeight','bold');
title('spectrogram of whitened data','FontWeight','bold');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');
hold off;


%% plot a basic power spectrum for each

L = length(signal);
fftAudioData = fft(signal);

P2 = abs(fftAudioData/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure(); 
set(gcf, 'Color', 'k'); 
plot(f, P1);
title('Power Spectrum of Raw Audio Data');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
xlim([0 120]);
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');

% Design a notch filter to remove 60 Hz noise
wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_60 = wo_60/5;      % Bandwidth
[b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
filteredAudio_1 = filter(b, a, signal_fragment);

L = length(filteredAudio_1);
fftAudioData = fft(filteredAudio_1);

P2 = abs(fftAudioData/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure(); 
set(gcf, 'Color', 'k'); 
plot(f, P1);
title('Power Spectrum of Bandpass Audio Data');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
xlim([0 120]);
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');

% FFT the signal and 
fftSignal = fft(filteredAudio_1);
powerSpectrum = abs(fftSignal).^2;
whitenedFFTSignal = fftSignal ./ sqrt(powerSpectrum + eps); % eps is added for numerical stability
whitenedSignal = real(ifft(whitenedFFTSignal));
whitenedSignal = whitenedSignal / max(abs(whitenedSignal));

L = length(whitenedSignal);
fftAudioData = fft(whitenedSignal);

P2 = abs(fftAudioData/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure(); 
set(gcf, 'Color', 'k'); 
plot(f, P1);
title('Power Spectrum of Whitened Audio Data');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
xlim([0 120]);
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');

%% Examples of clicks etc 

batdate = 221126;
mic = 1;
chunk = 3;
Fs = 192000;
    
%% Load in the audio data
audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');
chunk_filename = dir(strcat(audio_base,'/audioConCat_mic_',num2str(mic),'_segment_chunk_',num2str(chunk),'_break_*.mat'));
load(strcat(audio_base,chunk_filename.name));

seg_start_in_chunk = 30000000;
seg_end_in_chunk = 40000000;

signal = audioConCat(seg_start_in_chunk:seg_end_in_chunk);

% Feeder click
ssic = seg_start_in_chunk+18.45*Fs; % feeder click (chunk 3; batdate = 221126; mic = 1; seg_start_in_chunk = 30000000; seg_end_in_chunk = 40000000;
seic = seg_start_in_chunk+19.62*Fs; % feeder click ^^

signal_fragment = audioConCat(ssic:seic);

% Design a notch filter to remove 60 Hz noise
wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_60 = wo_60/5;      % Bandwidth
[b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
filteredAudio_1 = filter(b, a, signal_fragment);

L = length(filteredAudio_1);
fftAudioData = fft(filteredAudio_1);

P2 = abs(fftAudioData/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure(); 
set(gcf, 'Color', 'k'); 
plot(f, P1);
title('Power Spectrum of Feeder Click');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
xlim([0 400]);
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');


% Echolocations
ssic = seg_start_in_chunk+45.5*Fs; % echolocation train (chunk 3; batdate = 221126; mic = 1; seg_start_in_chunk = 30000000; seg_end_in_chunk = 40000000;
seic = seg_start_in_chunk+47.5*Fs; % echolocation train ^^

signal_fragment = audioConCat(ssic:seic);

% Design a notch filter to remove 60 Hz noise
wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_60 = wo_60/5;      % Bandwidth
[b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
filteredAudio_1 = filter(b, a, signal_fragment);

L = length(filteredAudio_1);
fftAudioData = fft(filteredAudio_1);

P2 = abs(fftAudioData/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure(); 
set(gcf, 'Color', 'k'); 
plot(f, P1);
title('Power Spectrum of Echo Train');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
xlim([0 400]);
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');

batdate = 221126;
mic = 1;
chunk = 5;
Fs = 192000;
    
%% Load in the audio data
audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');
chunk_filename = dir(strcat(audio_base,'/audioConCat_mic_',num2str(mic),'_segment_chunk_',num2str(chunk),'_break_*.mat'));
load(strcat(audio_base,chunk_filename.name));

seg_start_in_chunk = 30000000;
seg_end_in_chunk = 40000000;

signal = audioConCat(seg_start_in_chunk:seg_end_in_chunk);

% Human clicks
ssic = 159797612-(192000*3.6); % human click train (chunk 5; batdate = 221126; mic = 1)
seic = 159797612+192000*1; % human click train ^^

signal_fragment = audioConCat(ssic:seic);

% Design a notch filter to remove 60 Hz noise
wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_60 = wo_60/5;      % Bandwidth
[b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
filteredAudio_1 = filter(b, a, signal_fragment);

L = length(filteredAudio_1);
fftAudioData = fft(filteredAudio_1);

P2 = abs(fftAudioData/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure(); 
set(gcf, 'Color', 'k'); 
plot(f, P1);
title('Power Spectrum of Human Clicks');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
xlim([0 400]);
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');


% Human voice
ssic = 160871867; % madeleine "Ceasar" (chunk 5; batdate = 221126; mic = 1)
seic = 160871867+192000*1.5; % madeleine "Ceasar"

signal_fragment = audioConCat(ssic:seic);

% Design a notch filter to remove 60 Hz noise
wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_60 = wo_60/5;      % Bandwidth
[b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
filteredAudio_1 = filter(b, a, signal_fragment);

L = length(filteredAudio_1);
fftAudioData = fft(filteredAudio_1);

P2 = abs(fftAudioData/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure(); 
set(gcf, 'Color', 'k'); 
plot(f, P1);
title('Power Spectrum of Human Madeleine');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
xlim([0 400]);
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');






       