function [sorted_events_zerod,pk_loc_vec_zerod] = HumanBat_align_audio_to_ciholas(Bat,batdate,logger,alignment_,sorted_events,pk_loc_vec,first_ttl_sample,last_ttl_sample)

    % Sampling rates
    motu_Fs = 192000; % audio sampling rate
    ciholas_Fs = 120; % ciholas sampling rate
    
    exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/');
    
    % Load in the flight data 
    load(strcat(exp_data_path,'ciholas/aligned_bat_position_data_',num2str(Bat),'.mat'));
    
    % Load in the Extracted Behavior data
    load(strcat(exp_data_path,'ciholas/Extracted_Behavior_', num2str(batdate),'_',num2str(Bat),'.mat'));
    
    % Load in spike data
    load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
    
    if strcmp(alignment_,'ending')
        % When aligning using the LAST ttl, 
        % Add buffer to the start. The buffer is how much longer the ciholas_r vector is compared to the audio ttl vector
        gap_to_fill = (length(ciholas_r)/ciholas_Fs) - (last_ttl_sample/motu_Fs);
        % Scale this new gap by the audio sampling rate (motu samping rate)
        gap_to_fill_motu = gap_to_fill*motu_Fs;
        
        % Add this buffer to the ttl vector 
        pk_loc_vec_zerod = pk_loc_vec + gap_to_fill_motu;
        % Scale this ttl vector to the ciholas sample rate 
        pk_loc_vec_zerod_cih = round(pk_loc_vec_zerod / motu_Fs * ciholas_Fs);
    
        %% Load in the sorted_event times and add the same amount
        sorted_events_zerod = sorted_events + gap_to_fill_motu;
        % Sample the sorted_events_zerod to seconds, and then to samples 
        sorted_events_zerod_cih = round(sorted_events_zerod / motu_Fs * ciholas_Fs);
    
    elseif strcmp(alignment_,'start')
    
        % If you want to align to the START... subtract the first ttl from all values.
        
        pk_loc_vec_zerod = pk_loc_vec - first_ttl_sample+1;
        pk_loc_vec_zerod_cih = round(pk_loc_vec_zerod / motu_Fs * ciholas_Fs);
    
        %% Load in the sorted_event times and add the same amount
        sorted_events_zerod = sorted_events - first_ttl_sample + 1;
        % Sample the sorted_events_zerod to seconds, and then to samples 
        sorted_events_zerod_cih = round(sorted_events_zerod / motu_Fs * ciholas_Fs);       
    end
end