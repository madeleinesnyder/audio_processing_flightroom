%% Filter wav file

signal = data;
Fs = 192000;

% Design a notch filter to remove 60 Hz noise
wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
bw_60 = wo_60/2;      % Bandwidth
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

%% Write re-filtered wav

filename_audio = 'trial41_8.2286_echo1_bp.wav';
audioData = filteredAudio_2;
Fs = 192000;
audiowrite(filename_audio,audioData,Fs);

