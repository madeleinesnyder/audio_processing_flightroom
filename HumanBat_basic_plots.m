%% Make basic plots for a segment of unfiltered audio data 
function [] = HumanBat_basic_plots(signal,fs,sound_type,onefig)

    Fs = fs; motu_Fs = fs;
    
    %% Apply filters and repeat the process
    
    % % Design a notch filter to remove 60 Hz noise
    % wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
    % bw_60 = wo_60/5;      % Bandwidth
    % [b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
    % filteredAudio_1 = filter(b, a, signal);
    
    % Design a bandpass for what kind you want
    if strcmp(sound_type,'madeleine') | strcmp(sound_type,'kevin')
    
        % Design a notch filter to remove 60 Hz noise
        wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
        bw_60 = wo_60/5;      % Bandwidth
        [b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
        filteredAudio_1 = filter(b, a, signal);
    
        % bandpass from 10kHz to 40kHz
        fpass_lo = 40;
        fpass_hi = 10000;
        db_noise = 80;
        fband = 50;
    
         % Bandpass filter
        [B,A] = butter(2,[fpass_lo fpass_hi]/(motu_Fs/2));
        signal_bp = filtfilt(B,A,filteredAudio_1);
    
    elseif strcmp(sound_type,'feeder')
        % bandpass for feeder
        fpass_lo = 2000; 
        fpass_hi = 10000; 
        db_noise = 80;
        fband = 50;
            
        % Bandpass filter
        [B,A] = butter(2,[fpass_lo fpass_hi]/(motu_Fs/2));
        signal_bp = filtfilt(B,A,signal);
       
    elseif strcmp(sound_type,'echo')
        % bandpass from 20kHz to 50kHz

         % Design a notch filter to remove 60 Hz noise
        wo_60 = 60/(Fs/2);  % 60 Hz frequency in normalized frequency units
        bw_60 = wo_60/5;      % Bandwidth
        [b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
        filteredAudio_1 = filter(b, a, signal);

        fpass_lo = 40000; 
        fpass_hi = 60000; 
        db_noise = 90;
        fband = 100;

        fpass_lo = 10000;
        fpass_hi = 40000;
    
        % Bandpass filter
        [B,A] = butter(2,[fpass_lo fpass_hi]/(motu_Fs/2));
        signal_bp = filtfilt(B,A,filteredAudio_1);
    
    elseif strcmp(sound_type,'humanclick')
        % bandpass from 10kHz to 40kHz
        fpass_lo = 8000;
        fpass_hi = 30000;
        db_noise = 80;
        fband = 100;
    
         % Bandpass filter
        [B,A] = butter(2,[fpass_lo fpass_hi]/(motu_Fs/2));
        signal_bp = filtfilt(B,A,signal);
    end
    
    % Make spectrogram using Julie's code.
    if onefig == 0
        [to, fo, logB, pg, tError, fError] = HumanBat_Julie_Spec_Only_Bats(signal_bp, Fs, db_noise, fpass_hi, fband);  
    else
        [to, fo, logB, pg, tError, fError] = HumanBat_Julie_Spec_Only_Bats_onefig(signal_bp, Fs, db_noise, fpass_hi, fband);  
    end

    %% Get rid of 60Hz noise
    wo_60 = 60/(motu_Fs/2);  % 60 Hz frequency in normalized frequency units
    bw_60 = wo_60/5;      % Bandwidth
    [b, a] = iirnotch(wo_60, bw_60);  % IIR Notch filter design
    filteredAudio_1 = filter(b, a, signal);
    
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

end