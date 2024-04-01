function [all_audio_events] = bootleg_prune_echo_samples(all_audio_events)

    % time between samples min
    limit = 120;

    %% Function to eliminate all audio events that are too close together 
    diff_ae = [NaN,diff(all_audio_events)]*1000;
    sub_limit_idxs = find(diff_ae <= limit);
    
    % Eliminate all indexes that are duplicate echolocations
    all_audio_events(sub_limit_idxs) = [];

end
