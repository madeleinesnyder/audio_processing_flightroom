%% For a given session, find all the audio events and save them in a .mat file 

% Takes a reasonable amount of time.
% Run after HumanBat_batch_concatenate_audio.m 

% to look at: 260574388

% Change the bat date
batdate = 221126;

% Hyperparameters
Fs = 192000;
mic_num = 4;

%% Load in the pre-concatenated audio data
audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');
file_parts = dir(audio_base); 
file_parts = strsplit(file_parts(end-3).name,'_');
max_chunks = str2num(file_parts{end-2});

running_chunk_start = zeros(1,mic_num);
event_locs = cell(1,mic_num);
for mm = 1:mic_num
    for tt = 1:max_chunks

        % Find the timestamps of audio events (scaled_zerod_loc)
        % Find the length (in sampels) of the chunk of audio 
        viz = 0; % Set to 1 if you want to plot all the segments and examine
        seg_size = 5000000;
        [scaled_zerod_loc,len_chunk] = HumanBat_find_echo_times(batdate,mm,tt,seg_size,viz);

        AA = running_chunk_start(mm) + cell2mat(scaled_zerod_loc);
    
        % Add the length of the chunk of audio to the total amount of length in the chunks so far
        running_chunk_start(mm) = running_chunk_start(mm) + len_chunk;
        event_locs{mm} = [event_locs{mm},AA];
    end
end
% Save the timestamps cell array
save(strcat(audio_base,'event_timestamps.mat'),'event_locs');

%% Knit together the echolcoation clicks from all the microphones. Plot samples!
event_vec = [];
for mm = 1:mic_num
    event_vec = [event_vec,event_locs{mm}];
end

% Sort all audible events from all microphones chronologically and remove any duplicates. 
sorted_events = unique(sort(event_vec));
% Save the sorted, unique audio events. (Audio Sample Land)
save(strcat(audio_base,'sorted_event_timestamps.mat'),'sorted_events');

