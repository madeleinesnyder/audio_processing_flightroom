%% Script to downsample the audio in aprocessed folder to 16kHz for feeding into the NN

% Load in the segment of data
% Change the bat date
batdate = 220412;
Round = '2992814650';

% Hyperparameters
Fs = 192000;
mic_num = 4;

%% Load in the pre-concatenated audio data
audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/',Round,'/processed/',num2str(batdate),'/audio/audioConCat*chunk*');
file_parts = dir(audio_base); 
if batdate==221130
    max_chunks=7;
elseif batdate==221126
    max_chunks=8;
elseif batdate==221125
    max_chunks=8;
elseif batdate==221127
    max_chunks=8;
elseif batdate==221128
    max_chunks=8;
elseif batdate==221129
    max_chunks=7;
elseif batdate==221202
    max_chunks=7;
elseif batdate==221203
    max_chunks=6;
elseif batdate==220404
    max_chunks=7;
elseif batdate==220406
    max_chunks=8;
elseif batdate==220407
    max_chunks = 7;
elseif batdate==220408
    max_chunks = 8;
elseif batdate==220411 
    max_chunks = 6;
elseif batdate==220412
    max_chunks = 7;
elseif batdate==220413
    max_chunks = 0;
end

% Downsample it to 16kHz
originalSampleRate = 192000;
newSampleRate = 16000;
for i=1:max_chunks
    for jj = 1:mic_num
        audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/',Round,'/processed/',num2str(batdate),'/audio/audioConCat_mic_',num2str(jj),'_segment_chunk_',num2str(i),'*');
        file_parts = dir(audio_base); 
        name_to_load_done = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/',Round,'/processed/',num2str(batdate),'/audio/',file_parts(end).name(1:41),'_16kHz.wav');
        if exist(name_to_load_done)
            continue;
        else
            name_to_load = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/',Round,'/processed/',num2str(batdate),'/audio/',file_parts(1).name);
            load(name_to_load);
            
            % Calculate the greatest common divisor to find the resampling factors
            gcdSampleRate = gcd(originalSampleRate, newSampleRate);
            p = newSampleRate / gcdSampleRate;
            q = originalSampleRate / gcdSampleRate;
            
            % Resample the audio to the new sample rate
            downsampledAudio = resample(audioConCat, p, q);
            
            % Write the downsampled audio to a new file
            save_name = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/',Round,'/processed/',num2str(batdate),'/audio/',file_parts(1).name(1:end-4),'_16kHz.wav');
            audiowrite(save_name, downsampledAudio, newSampleRate);
        end
    end
end

close all; clearvars -except max_chunks mic_num batdate newSampleRate Round

% Load in all 16Hz samples and concatentate into one
% Downsample it to 16kHz
for jj = 1:mic_num
    audio_part=cell(1,max_chunks);
    for i=1:max_chunks
        file_parts = dir(strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/',Round,'/processed/',num2str(batdate),'/audio/audioConCat_mic_',num2str(jj),'_segment_chunk_',num2str(i),'*_16kHz.wav'));
        name_to_load = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/',Round,'/processed/',num2str(batdate),'/audio/',file_parts.name);
        audio_part{i} = audioread(name_to_load);
    end
    audio_full = cell2mat(audio_part');
    save_name =strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/',Round,'/processed/',num2str(batdate),'/audio/',file_parts.name(1:18),'16kHz.wav');
    audiowrite(save_name, audio_full, newSampleRate);
end

% Load in all mics and concat to one 
all_mics=cell(1,mic_num);
for jj = 1:mic_num
    load_name =strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/',Round,'/processed/',num2str(batdate),'/audio/audioConCat_mic_',num2str(jj),'_16kHz.wav');
    all_mics{jj} = audioread(load_name);
end
allmics_full = cell2mat(all_mics');
save_name = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/',Round,'/processed/',num2str(batdate),'/audio/audioConCat_all_mics_16kHz.wav');
audiowrite(save_name, allmics_full, newSampleRate);
