%% Make basic plots for a segment of unfiltered audio data 
function [] = HumanBat_basic_plots(signal,fs,sound_type)

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
        fpass_hi = 2000;
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
        fpass_lo = 40000; 
        fpass_hi = 60000; 
        db_noise = 80;
        fband = 100;
    
        % Bandpass filter
        [B,A] = butter(2,[fpass_lo fpass_hi]/(motu_Fs/2));
        signal_bp = filtfilt(B,A,signal);
    
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
    [to, fo, logB, pg, tError, fError] = HumanBat_Julie_Spec_Only_Bats(signal_bp, Fs, db_noise, fpass_hi, fband);  
end