% Concatenate the raw audio into 100-"audio-trial" segemnts and save in
% that date's processed directory. 

%% Audio extraction pipeline 

% All scripts for the audio extraction pipeline are in:
% /home/madeleine/Desktop/HumanSniffBat/HumanBat/audio/batch_process

% 1. Batch extract the raw audio. HumanBat_batch_concatenate_audio.m (this script!).
% % Make sure all [date]_audio_trial_micX_batch_0.mat files are present in the directory before running
% subsequent steps!

% 2. Find the samples at which audio events occur. HumanBat_find_echo_times_audioclip.m.
% The samples are the sample in the context of the whole session (batdate), concatenation accounted for.

% 3. Find the samples at which TTL pulses occur. HumanBat_get_ttl_peaks.m.
% Extra samples at the beginning and end are accounted for. 

% 4. Plot X seconds around a desired audio event and plot aligned to
% flight. HumanBat_visualize_audio_event_aligned_to_flight.m

% 5. Plot 50ms around all audio events, human-label them, segment and store
% them for BioSound processing. HumanBat_shortclip_label_and_save.m

batdate_list = [221203];

% For every session in the list...
for i=1:length(batdate_list)
    batdate = batdate_list(i);
    disp(strcat("Copying audio for ",num2str(batdate)));

    chdir(strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate)));
    if ~exist("audio")
        system('sudo mkdir audio');
    end
    
    % Copy raw audio files to "processed" directory 
    source_filepath = strcat(strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/raw/',num2str(batdate),'/audio/*'));
    [file_parts] = dir(source_filepath); source_filepath = strcat(source_filepath(1:end-2),'/',file_parts(end).name);
    system(strcat("cp -r "," ",source_filepath," ",strcat(pwd,"/audio")));
    
    % Concatenate the raw audio in the processed directory into 100-trial segments
    disp(strcat("Concatenating audio for ",num2str(batdate)));
    files = dir("audio");
    chdir(strcat("audio/"));
    HumanBat_audio_concat(pwd);
    
    % Remove raw files from processed folder 
    disp(strcat("Removing audio trial files from processed dir for ",num2str(batdate)));
    raw_files = dir(strcat(pwd,'/_',num2str(batdate),'*'));
    name_ = strcat('_',num2str(batdate),'_audio_trial_*');
    arg2 = strcat("find -type f -name '",name_,"' -delete");
    system(arg2);
    
end
    

