% Ask which bat an echolocation belongs to 
function [useable_echo_timestamps,who_echod,probably_not_echos,mic_list] = bootleg_who_echolocated(bat1_ciholas_r,bat2_ciholas_r,mic_audio_events,valid_flight_samples)

    % For every audio event, extract the events that occur during these
    % flights
    all_audio_events = horzcat(cell2mat(mic_audio_events));
    useable_echo_timestamps = [];
    mic_list = {}; % Contains which mic this echolocation was closest to 
    echo_dist_thresh = 2500;

    % For each echolocation, assign a bat based on how close that bat was
    % to a given microphone 
    mic_positions = bootleg_hardcode_mic_positions();

    % For each echolocation, determine which bat was closest to the
    % microphone that picked it up 
    who_echod = []; who_echod_distances = [];
    for i=1:length(all_audio_events)
        if ismember(round(all_audio_events(i)*120),valid_flight_samples)
            useable_echo_timestamps = [useable_echo_timestamps,round(all_audio_events(i)*120)];
            which_mics = [];
            for j=1:length(mic_audio_events)
                if ismember(all_audio_events(i),mic_audio_events{j})
                    which_mics = [which_mics,j];
                end
            end
            if length(which_mics) > 1
                who_echod_distances_temp = nan(length(which_mics),2);
                disp("echolocation evenly detected between two microphones.")
                for m=1:length(which_mics)
                    mic_pos = mic_positions(which_mics(m),:);
                    bat1_from_mic = pdist2(mic_pos,bat1_ciholas_r(round(all_audio_events(i)*120),:));
                    bat2_from_mic = pdist2(mic_pos,bat2_ciholas_r(round(all_audio_events(i)*120),:));
                    who_echod_distances_temp(m,:) = [bat1_from_mic,bat2_from_mic];
                end
                if sum(who_echod_distances_temp(1,:) < who_echod_distances_temp(2,:)) == 2
                    echo_bat = 1;
                elseif sum(who_echod_distances_temp(1,:) < who_echod_distances_temp(2,:)) == 0
                    echo_bat = 2;
                else
                    echo_bat = NaN;
                end
                who_echod = [who_echod,echo_bat];
                find(sum(who_echod_distances_temp,2) == min(sum(who_echod_distances_temp,2)))
                who_echod_distances = [who_echod_distances;who_echod_distances_temp(find(sum(who_echod_distances_temp,2) == min(sum(who_echod_distances_temp,2))),:)];
            elseif isempty(which_mics)
                continue
            elseif length(which_mics) == 1
                % Which bat was closest to the mic at this time?
                mic_pos = mic_positions(which_mics,:);
                bat1_from_mic = pdist2(mic_pos,bat1_ciholas_r(round(all_audio_events(i)*120),:));
                bat2_from_mic = pdist2(mic_pos,bat2_ciholas_r(round(all_audio_events(i)*120),:));
                if bat1_from_mic < bat2_from_mic
                    echo_bat = 1;
                else
                    echo_bat = 2;
                end
                who_echod = [who_echod,echo_bat];
                who_echod_distances = [who_echod_distances;[bat1_from_mic,bat2_from_mic]];
            end
            mic_list{end+1} = which_mics;
        end
    end

    probably_not_echos = all(who_echod_distances>echo_dist_thresh,2);
 
end
