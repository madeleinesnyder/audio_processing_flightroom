%% For a given session, find all the audio events around a specific part of the data 
 %[ debugging file ]%

% Change the bat date
batdate = 221126;

% Hyperparameters
Fs = 192000;
mic_num = 1;

%% Load in the pre-concatenated audio data
audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');
file_parts = dir(audio_base); 
file_parts = strsplit(file_parts(end-3).name,'_');
max_chunks = str2num(file_parts{end-2});

running_chunk_start = zeros(1,mic_num);
event_locs = cell(1,mic_num);

% WHERE DO YOU WANT TO FIND AUDIO EVENTS AROUND?
desired_event = 771504467;
desired_event = 415439364

for mm = 1:mic_num
    for tt = 1:max_chunks

        % Find the timestamps of audio events (scaled_zerod_loc)
        % Find the length (in sampels) of the chunk of audio 
        viz = 0; % Set to 1 if you want to plot all the segments and examine
        seg_size = 5000000;
        [len_chunk] = HumanBat_find_len_chunk(batdate,mm,tt,seg_size,viz);
        running_chunk_start(mm) = running_chunk_start(mm) + len_chunk;

        if (desired_event < running_chunk_start(mm)) & (desired_event > (running_chunk_start(mm) - len_chunk))
            amt_into = desired_event - (running_chunk_start(mm)-len_chunk);
            [scaled_zerod_loc,len_chunk] = HumanBat_find_echo_times(batdate,mm,tt,seg_size,viz);
            AA = running_chunk_start(mm) + cell2mat(scaled_zerod_loc);
            event_locs{mm} = [event_locs{mm},AA];
        end
    end
end