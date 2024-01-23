%% Plot AA

batdate = 221126;
seg_start = 1;
seg_end=5000000;
mic_num=4;
mm=1
motu_Fs = 192000;

audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');

whitenedSignal_mm = {};
%for mm = 1:mic_num 
    
    audio_file = dir(strcat(audio_base,'/audioConCat_mic_',num2str(mm),'_segment_chunk_1_break_*.mat'));
    load(strcat(audio_base,audio_file.name));

    % Apply segment bounds to desired segment to plot 
    signal = audioConCat(seg_start:seg_end);
    
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

%end

whitenedSignal_allmics = cell2mat(whitenedSignal_mm);
sum_whitenedSignal_allmics = sum(whitenedSignal_allmics,2);
    
% above_zero_events = (sorted_events-event_vec_start);
% above_zero_events = above_zero_events(above_zero_events>0);
% above_zero_events = above_zero_events(above_zero_events<length(sum_whitenedSignal_allmics));
% above_zero_events_aligned_to_start_of_session = above_zero_events+event_vec_start;

%% Align the above_zero_events to the flight data 
above_zero_events = (event_locs-event_vec_start);
%above_zero_events = above_zero_events(above_zero_events>=0);
%above_zero_events(1) = 1;
%above_zero_events = above_zero_events(above_zero_events<length(sum_whitenedSignal_allmics));
%above_zero_events_aligned_to_start_of_session = above_zero_events+event_vec_start;
%[above_zero_events_aligned,ttl_events_aligned] = HumanBat_align_audio_to_ciholas(Bat,batdate,logger,alignment_,above_zero_events_aligned_to_start_of_session,pk_loc_vec,first_ttl_sample,last_ttl_sample);
load('scaled_loc.mat');
AA_less = scaled_loc(scaled_loc<length(signal));


figure(); 
subplot(4,1,1); hold on; 
title("All mic echolocations plotted on filtered audio",'FontWeight','bold','Color','w');
%plot(sum_whitenedSignal_allmics);
    plot(whitenedSignal_mm{mm});
ee = NaN(1,length(sum_whitenedSignal_allmics));
ee(AA_less) = 1;
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
%[S,F,T,P]=spectrogram(sum_whitenedSignal_allmics,w,overlap,nfft,motu_Fs);
%[S2]=spectrogram(sum_whitenedSignal_allmics,dw,overlap,nfft,motu_Fs);
[S,F,T,P]=spectrogram(whitenedSignal_mm{1},w,overlap,nfft,motu_Fs);
[S2]=spectrogram(whitenedSignal_mm{1},dw,overlap,nfft,motu_Fs);
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

spectrogram_ee_idxs = floor(AA_less/(length(sum_whitenedSignal_allmics)/length(T)));
spectrogram_ee = NaN(1,length(T));
if spectrogram_ee_idxs(1) == 0; spectrogram_ee_idxs(1) = 1; end;
spectrogram_ee(spectrogram_ee_idxs) = 60000;

% Plot the spectrogram with the peaks highlighted
subplot(4,1,2); hold on; 
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gcf, 'Color', 'k'); 
colormap(hot)
imagesc(T,F,IMAGE_MODS); %set(gca,'YDir','normal');
scatter(T, spectrogram_ee',5,'y'); % Mark the peaks
ylabel('Frequency kHz','FontWeight','bold');
xlabel('time (s)','FontWeight','bold');
title('Spectrogram of whitened data','FontWeight','bold','Color','w');
axis tight;
hold off;
