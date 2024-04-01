% Turn bootleg timestamps into .wav files

% load in session timestamps 

function [] = bootleg_make_wavs(Bat,batdate)

    load(strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/saved_echolocation_timestamps.mat'));
    save_name_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/saved_event_wav_');

    for i=1:4
        % Load in audio file
        mic = i;
        [audio_data,Fs_read_in] = audioread(strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat_KQ/1464314684/processed/',num2str(batdate),'/audio/Flight_Mic',num2str(mic),'_',num2str(batdate),'_.wav'));
        seg_len = 0.25*Fs_read_in;

        % Make wavs 
        for j=1:length(mic_audio_events{i})
            % Load in audio file and make the wavs
    
            y = audio_data(round(mic_audio_events{i}(j)*Fs_read_in)-seg_len:round(mic_audio_events{i}(j)*Fs_read_in)+seg_len);
            filename = strcat(save_name_base,'mic_',num2str(i),'_',num2str(mic_audio_events{i}(j)),'.wav');
    
            % Write the file
            audiowrite(filename, y, Fs_read_in);
        end
    end
    disp("Successfully wrote all .wav files!")
end


