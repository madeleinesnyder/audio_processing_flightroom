%% Load in the ds audio data

mic_audio_events = {};
for mm=1:4
    audio_events_s = [];
    mic = mm;
    Round = '1464314684';
    batdate = 221130;
    [audio_data,motu_Fs] = audioread(strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat_KQ/',Round,'/processed/',num2str(batdate),'/audio/Flight_Mic',num2str(mic),'_',num2str(batdate),'_.wav'));

    % chop out segment of data we want to look at
    num_splits = 48;
    for j=1:num_splits

        % Load in audio
        [audio_data,motu_Fs] = audioread(strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat_KQ/',Round,'/processed/',num2str(batdate),'/audio/Flight_Mic',num2str(mic),'_',num2str(batdate),'_.wav'));
        start_of_j = (j-1)*round(length(audio_data)/num_splits)+1;
        if j==num_splits
            seg = audio_data(start_of_j:end);
        else
            seg = audio_data(start_of_j:start_of_j+round(length(audio_data)/num_splits));
        end
        clear audio_data;
    
        %% Get rid of 60Hz noise
        wo_60 = 60/(motu_Fs/2);  % 60 Hz frequency in normalized frequency units
        bw_60 = wo_60/5;      % Bandwidth
        [b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
        filteredAudio_1 = filter(b, a, seg);
        
        %% FFT the signal and whiten
        fftSignal = fft(filteredAudio_1);
        powerSpectrum = abs(fftSignal).^2;
        
        whitenedFFTSignal = fftSignal ./ sqrt(powerSpectrum + eps); % eps is added for numerical stability
        whitenedSignal = real(ifft(whitenedFFTSignal));
        
        whitenedSignal_mm= whitenedSignal / max(abs(whitenedSignal));
        whitenedSignal_mm(1:10) = 0; whitenedSignal_mm(end-10:end) = 0;
    
        %% Hi pass filter
        cutoffFrequency = 4000; % Cutoff frequency in Hz
        filterOrder = 2; % Order of the filter
        
        % Step 2: Design the filter
        % Design a Butterworth lowpass filter
        d = designfilt('highpassiir', ...
                       'FilterOrder', filterOrder, ...
                       'HalfPowerFrequency', cutoffFrequency, ...
                       'SampleRate', motu_Fs);
        
        % Step 3: Apply the filter
        seg_y = filter(d, whitenedSignal_mm);
    
        %% Make figure
    
%         figure(); 
%         subplot(4,1,1); hold on; 
%         title("All mic echolocations plotted on filtered audio",'FontWeight','bold','Color','w');
%         plot(seg_y);
%         %ee_echo = NaN(1,length(whitenedSignal_mm)); ee_feeder = NaN(1,length(whitenedSignal_mm)); ee_humanclick = NaN(1,length(whitenedSignal_mm));
%         %ee_echo(echo_events) = 1; ee_feeder(feeder_events) = 1; ee_humanclick(humanclick_events) = 1;
%         %scatter([1:length(sum_whitenedSignal_allmics)],ee_echo,5,'filled','g');
%         %scatter([1:length(sum_whitenedSignal_allmics)],ee_humanclick,5,'filled','y');
%         %scatter([1:length(sum_whitenedSignal_allmics)],ee_feeder,5,'filled','w');
%         set(gca, 'Color', 'k');
%         set(gca, 'XColor', 'w');
%         set(gca, 'YColor', 'w');
%         set(gcf, 'Color', 'k'); 
%         xlabel("time (samples)",'FontWeight','bold');
%         ylabel("amplitude",'FontWeight','bold');
%         axis tight;
%         hold off;
%         
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
        [S,F,T,P]=spectrogram(seg_y,w,overlap,nfft,motu_Fs);
        [S2]=spectrogram(seg_y,dw,overlap,nfft,motu_Fs);
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
        
        %spectrogram_ee_echo = NaN(length(T),1); spectrogram_ee_feeder = NaN(length(T),1); spectrogram_ee_hc = NaN(length(T),1);
        %spectrogram_ee_echo_idxs = floor(echo_events/(length(sum_whitenedSignal_allmics)/length(T))); spectrogram_ee_feeder_idxs = floor(feeder_events/(length(sum_whitenedSignal_allmics)/length(T))); spectrogram_ee_hc_idxs = floor(humanclick_events/(length(sum_whitenedSignal_allmics)/length(T)));
        %if ~isempty(echo_events)
        %    if spectrogram_ee_echo_idxs(1) == 0; spectrogram_ee_echo_idxs(1) = 1; end;
        %end
        %if ~isempty(feeder_events)
        %    if spectrogram_ee_feeder_idxs(1) == 0; spectrogram_ee_feeder_idxs(1) = 1; end;
        %end
        %if ~isempty(humanclick_events)
        %    if spectrogram_ee_hc_idxs(1) == 0; spectrogram_ee_hc_idxs(1) = 1; end;
        %end
        %spectrogram_ee_echo(spectrogram_ee_echo_idxs) = 60000; spectrogram_ee_feeder(spectrogram_ee_feeder_idxs) = 60000; spectrogram_ee_hc(spectrogram_ee_hc_idxs) = 60000;
        
%         % Plot the spectrogram with the peaks highlighted
%         subplot(4,1,2); hold on; 
%         colormap(hot)
%         imagesc(T,F,IMAGE_MODS); set(gca,'YDir','normal');
%         axis tight;
%         %plot(T, spectrogram_ee_echo','g*','MarkerSize',5); % Mark the peaks
%         %plot(T, spectrogram_ee_feeder','w*','MarkerSize',5); % Mark the peaks
%         %plot(T, spectrogram_ee_hc','y*','MarkerSize',5); % Mark the peaks
%         ylabel('Frequency kHz','FontWeight','bold');
%         xlabel('time (s)','FontWeight','bold');
%         title('Spectrogram of whitened data','FontWeight','bold','Color','w');
%         set(gca, 'Color', 'k');
%         set(gca, 'XColor', 'w');
%         set(gca, 'YColor', 'w');
%         set(gcf, 'Color', 'k'); 
%         hold off;
        
        sum_IMAGE_MODS = sum(IMAGE_MODS,1);
        
        %% Findpeaks 
        
        min_peak_dist = 0.02/96000*size(IMAGE_MODS,2)*10000;
        min_peak_height = 10;
        [putative_peaks,putative_idxs] = findpeaks((sum_IMAGE_MODS-min(sum_IMAGE_MODS)),'MinPeakDistance',min_peak_dist,'MinPeakHeight',min_peak_height);
        peak_vec = NaN(1,length(sum_IMAGE_MODS));
        peak_vec(putative_idxs) = 20;
        
        figure();
        subplot(2,1,1);
        colormap(hot)
        plot((sum_IMAGE_MODS-min(sum_IMAGE_MODS)));
        set(gca,'YDir','normal'); hold on;
        axis tight;
        scatter([1:length(sum_IMAGE_MODS)], peak_vec',50,'y*'); % Mark the peaks
        title('Spectrogram of whitened data','FontWeight','bold','Color','w');
        set(gca, 'Color', 'k');
        set(gca, 'XColor', 'w');
        set(gca, 'YColor', 'w');
        set(gcf, 'Color', 'k'); 
        hold off;
        
        peak_vec = NaN(1,length(sum_IMAGE_MODS));
        peak_vec(putative_idxs) = 40000;
        subplot(2,1,2);
        colormap(hot)
        imagesc(T,F,IMAGE_MODS); set(gca,'YDir','normal'); hold on;
        axis tight;
        scatter(T, peak_vec',50,'y*'); % Mark the peaks
        ylabel('Frequency kHz','FontWeight','bold');
        xlabel('time (s)','FontWeight','bold');
        title('Spectrogram of whitened data','FontWeight','bold','Color','w');
        set(gca, 'Color', 'k');
        set(gca, 'XColor', 'w');
        set(gca, 'YColor', 'w');
        set(gcf, 'Color', 'k'); 
        hold off;

        total_s = length(seg_y)/motu_Fs;
        scale_factor = round(length(sum_IMAGE_MODS)/(length(seg_y)/motu_Fs));

        % Scale the indexes of audio events
        scaled_events_s = putative_idxs/scale_factor;
        audio_events_s = [audio_events_s,scaled_events_s+start_of_j/motu_Fs];

        clearvars -except audio_events_s motu_Fs mm j Round batdate mic_audio_events num_splits mic
        close all;
    end
    
    mic_audio_events{end+1} = audio_events_s;
end

save(strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/saved_echolocation_timestamps.mat'),'mic_audio_events');

    
    
    
    

