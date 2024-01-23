function [len_chunk] = HumanBat_find_len_chunk(batdate,mic,chunk,segment_size,visualize)
    
    % Define audio sample rate 
    Fs = 192000;
    
    %% Load in the audio data
    audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');
    chunk_filename = dir(strcat(audio_base,'/audioConCat_mic_',num2str(mic),'_segment_chunk_',num2str(chunk),'_break_*.mat'));
    load(strcat(audio_base,chunk_filename.name));
    
    len_chunk = length(audioConCat);
end