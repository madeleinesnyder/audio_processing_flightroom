%% Align the audio events to the flights/ephys 

% Run after HumanBat_find_echo_times_audioclip.m

batdate = 221126;
Fs = 192000;
mic_num = 4;

%% Load in the ttl data to find the first TTL, last TTL, and TTL timestamps in milliseconds
audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/');

pk_loc_vec = [];
running_chunk = 0; 
for i=1:8
    chunk_filename = dir(strcat(audio_base,'ttlConCat_segment_chunk_',num2str(i),'_break_*.mat'));
    load(strcat(audio_base,chunk_filename.name));
    
    min_pk_dist = 500000;
    min_pk_ht = 0.35;
    
    [pk_vals,pk_locs] = findpeaks(ttlConCat,'MinPeakDistance',min_pk_dist,'MinPeakHeight',min_pk_ht);
    pk_loc_vec = [pk_loc_vec;running_chunk+pk_locs];
    
    % This is the value you will subtract from all vectors in audio space to
    % zero it to the first TTL.
    if i==1
        first_ttl_sample = running_chunk+pk_locs(1);
    end
    if i==8
        last_ttl_sample = running_chunk+pk_locs(end);
    end

    running_chunk = running_chunk + length(ttlConCat);
end

save(strcat(audio_base,'ttl_locs.mat'),'pk_loc_vec');
save(strcat(audio_base,'ttl_first_sample.mat'),'first_ttl_sample');
save(strcat(audio_base,'ttl_last_sample.mat'),'last_ttl_sample');

%% Check to see if there is a TTL in the last-last file

raw_audio_base_dir = dir(strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/raw/',num2str(batdate),'/audio/*'));
raw_audio_base = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/raw/',num2str(batdate),'/audio/',raw_audio_base_dir(end).name);
a_files = dir(raw_audio_base);
for i=1:20
    load(strcat(raw_audio_base,'/_',num2str(batdate),'_audio_trial_',num2str(size(a_files,1)-3-i),'.mat'));
    [rnd_vals,rnd_pks] = findpeaks(recbuf(:,end),'MinPeakProminence',0.3);
    if ~isempty(rnd_pks)
        disp("add X TTLs to the final TTL");
        for j=1:length(rnd_pks)
            pk_loc_vec = [pk_loc_vec;pk_loc_vec(end)+(pk_loc_vec(end)-pk_loc_vec(end-1))];
        end
        last_ttl_sample = pk_loc_vec(end);
        break;
    else
        disp("All good");
    end
end

% Re-save the values for the ttl locations and the last sample.
save(strcat(audio_base,'ttl_locs.mat'),'pk_loc_vec');
save(strcat(audio_base,'ttl_first_sample.mat'),'first_ttl_sample');
save(strcat(audio_base,'ttl_last_sample.mat'),'last_ttl_sample');






