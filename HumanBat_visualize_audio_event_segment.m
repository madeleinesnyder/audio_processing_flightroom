%% Plot the X ms around a given "audio event"

batdate = 221125;
Fs = 192000;
mic_num = 4;

%% Load in the audio data
audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');
load(strcat(audio_base,'event_timestamps.mat'),'event_locs');
load(strcat(audio_base,'sorted_event_timestamps.mat'),'sorted_events');

% Plot an example
chunk = 5; 
mic = 3;
seg_start = 135000000; seg_end = 140000000;

% Load in every chunk and get the final accumulated run start sample, then
% add the start and end seg to it.
run_chunk = 0;
for i=1:chunk-1

    audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');
    load(strcat(audio_base,'/audioConCat_mic_',num2str(mic),'_segment_chunk_',num2str(i),'_break_0.mat'));
    
    len_chunk = length(audioConCat);
    run_chunk = run_chunk+len_chunk;
end

if chunk == 1
    event_vec_start = seg_start;
    event_vec_end = seg_end;
else
    event_vec_start = run_chunk+seg_start;
    event_vec_end = run_chunk+seg_end;
end

%% Get filtered data for segment from all mics 
whitenedSignal_mm = {};
for mm = 1:mic_num 
    
    load(strcat(audio_base,'/audioConCat_mic_',num2str(mm),'_segment_chunk_',num2str(chunk),'_break_0.mat'));

    % Apply segment bounds to desired segment to plot 
    signal = audioConCat(seg_start:seg_end);
    
    %% Remove 60 and 90 Hz noise from segment
    
    % Design a notch filter to remove 60 Hz noise
    wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
    bw_60 = wo_60/5;      % Bandwidth
    [b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
    filteredAudio_1 = filter(b, a, signal);
        
    % Notch filter for 90 Hz
    wo_90 = 90/(Fs/2);  % 60 Hz frequency in normalized frequency units
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
    
above_zero_events = (sorted_events-event_vec_start);
above_zero_events = above_zero_events(above_zero_events>0);
above_zero_events = above_zero_events(above_zero_events<length(sum_whitenedSignal_allmics));

%% Plot the signal segment and the echolocations from each microphone

figure(); 
subplot(2,1,1); hold on; 
title("All mic echolocations plotted on filtered audio",'FontWeight','bold');
for mm=1:mic_num
    plot(whitenedSignal_mm{mm});
end
ee = NaN(1,length(sum_whitenedSignal_allmics));
ee(above_zero_events) = 1;
scatter([1:length(sum_whitenedSignal_allmics)],ee,10,'r');
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

t=-N/2+1:N/2;
sigma=(tscale/1e3)*Fs;
w = exp(-(t/sigma).^2);
dw = -2*w.*(t/(sigma^2));
    
% Calculate spectrogram
[S,F,T,P]=spectrogram(sum_whitenedSignal_allmics,w,overlap,nfft,Fs);
[S2]=spectrogram(sum_whitenedSignal_allmics,dw,overlap,nfft,Fs);
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

spectrogram_ee_idxs = round(above_zero_events/(length(sum_whitenedSignal_allmics)/length(T)));
spectrogram_ee = NaN(1,length(T));
spectrogram_ee(spectrogram_ee_idxs) = 60000;

% Plot the spectrogram with the peaks highlighted
subplot(2,1,2); hold on; 
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gcf, 'Color', 'k'); 
colormap(hot)
imagesc(T,F,IMAGE_MODS); %set(gca,'YDir','normal');
scatter(T, spectrogram_ee',5,'y'); % Mark the peaks
ylabel('Frequency kHz','FontWeight','bold');
xlabel('time (s)','FontWeight','bold');
title('Spectrogram of whitened data','FontWeight','bold');
axis tight;
hold off;


    








