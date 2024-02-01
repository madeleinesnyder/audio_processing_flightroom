%% Load in a chunk of audio and save the timestamps of echolocations to a file
function [scaled_zerod_loc,len_chunk] = HumanBat_find_echo_times(batdate,mic,chunk,segment_size,visualize)
    
    % Define audio sample rate 
    Fs = 192000;
    
    %% Load in the audio data
    audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');
    chunk_filename = dir(strcat(audio_base,'/audioConCat_mic_',num2str(mic),'_segment_chunk_',num2str(chunk),'_break_*.mat'));
    load(strcat(audio_base,chunk_filename.name));
    
    len_chunk = length(audioConCat);
    
    %% FFT visualization and echolocation ID using small segments.
    
    addpath(genpath('/home/madeleine/Desktop/HumanSniffBat/shared_utils/FinchScope'));
    segment_num = ceil(len_chunk/segment_size);
    
    % Store the locations of the echolocations
    scaled_zerod_loc = {};
    for ss = 1:segment_num
        
        echo_segment_start = (ss-1)*segment_size+1;
        echo_segment_end = (ss-1)*segment_size+segment_size;

        if ss == segment_num
            signal = audioConCat(echo_segment_start:end); % Sine wave
        else
            signal = audioConCat(echo_segment_start:echo_segment_end); % Sine wave
        end

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
        
        %% Whiten the signal whole chunk (flattening the power spectrum)
        
        % FFT the signal and 
        fftSignal = fft(filteredAudio_2);
        powerSpectrum = abs(fftSignal).^2;
        
        whitenedFFTSignal = fftSignal ./ sqrt(powerSpectrum + eps); % eps is added for numerical stability
        whitenedSignal = real(ifft(whitenedFFTSignal));
        
        whitenedSignal = whitenedSignal / max(abs(whitenedSignal));
           
        %t = 0:1/Fs:1; % Time vector (1 second)
        
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
        
        % S: Short-time Fourier Transform
        % F: Frequency vector
        % T: Time vector
        % P: Power spectral density
    
        %% Get the peaks of the summed frequency bands
        
        IMAGE_MODS_SUM = sum(IMAGE_MODS);
        threshold = 5;
        peak_dist_threshold_seconds = 0.0012*Fs;
        p_height = 0.05;
        [idx, loc] = findpeaks(abs(IMAGE_MODS_SUM), 'MinPeakProminence', threshold,'MinPeakDistance',peak_dist_threshold_seconds,'MinPeakHeight',p_height);
        
        % Eliminate all found locs if the unfiltered signal is just
        % oscillitory noise
        if max(filteredAudio_2) < 0.03
            loc = [];
        end

        % Turn locations into a binary vector and scale the 1's to high frequency
        % for plotting purposes.
        %loc_signaltime = loc*(length(signal_no_echo_bp)/length(IMAGE_MODS_SUM))
        pf_binary = NaN(length(T),1); pf_binary(loc) = 1;
        pf_binary(find(pf_binary==1)) = 63000;
    
        % Save the times of the echolocation peaks in the long vector. 
        scaled_loc = round(loc*(length(whitenedSignal)/length(T)));     % Scale the locations of the echolocations according to the samples in the real data vector
        
        scaled_loc_binary = NaN(length(whitenedSignal),1); scaled_loc_binary(scaled_loc) = 0.3;
        
        if visualize == 1
            % Plot the spectrogram with the peaks highlighted
            figure(); 
            subplot(2,1,1); hold on; 
            colormap(hot)
            imagesc(T,F,IMAGE_MODS); set(gca,'YDir','normal');
            axis tight;
            plot(T, pf_binary, 'y*', 'MarkerSize', 1); % Mark the peaks
            ylabel('frequency kHz','FontWeight','bold');
            xlabel('time (s)','FontWeight','bold');
            title('spectrogram of whitened data','FontWeight','bold');
            set(gca, 'Color', 'k');
            set(gca, 'XColor', 'w');
            set(gca, 'YColor', 'w');
            set(gca, 'GridColor', 'w');
            hold off;
            
            subplot(2,1,2); 
            set(gcf, 'Color', 'k'); hold on;
            title("Whitened signal with echolocation peaks",'FontWeight','bold');
            xlabel("time (samples) at 192000Hz",'FontWeight','bold');
            ylabel("amplitude",'FontWeight','bold');
            colors = [linspace(0, 1, length(whitenedSignal))', linspace(0, 0.5, length(whitenedSignal))', linspace(1, 0, length(whitenedSignal))'];
            set(gca, 'Color', 'k');
            set(gca, 'XColor', 'w');
            set(gca, 'YColor', 'w');
            set(gca, 'GridColor', 'w');
            plot(whitenedSignal); 
            scatter(1:length(whitenedSignal),scaled_loc_binary,5,'r','filled');
            hold off;
    
            figure(); 
            set(gcf, 'Color', 'k'); hold on;
            title("Whitened signal with echolocation peaks",'FontWeight','bold');
            xlabel("time (samples) at 192000Hz",'FontWeight','bold');
            ylabel("amplitude",'FontWeight','bold');
            colors = [linspace(0, 1, length(whitenedSignal))', linspace(0, 0.5, length(whitenedSignal))', linspace(1, 0, length(whitenedSignal))'];
            set(gca, 'Color', 'k');
            set(gca, 'XColor', 'w');
            set(gca, 'YColor', 'w');
            set(gca, 'GridColor', 'w');
            axis off; colorbar;
            y_norm = (whitenedSignal - min(whitenedSignal)) / (max(whitenedSignal) - min(whitenedSignal));
            scatter(1:length(whitenedSignal),whitenedSignal,3,y_norm,'filled'); 
            scatter(1:length(whitenedSignal),scaled_loc_binary,10,'red','filled');
            hold off;

            userInput = input('Press n for next segment and q to quit: ','s');

            if strcmp(userInput,'n')
                disp("to next segment");
            elseif strcmp(userInput,'q')
                break;
            end
        end
        
        % Save the locs of the echolocations relative to the first sample in
        % the audiofile
        scaled_zerod_loc{ss} = scaled_loc + echo_segment_start-1;

       % if ss == 5
       %     disp("o")
       % end
        
        clear scaled_loc loc;
    end
end
    
    


