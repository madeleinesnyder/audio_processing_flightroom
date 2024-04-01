function [] = bootleg_plot_echos_on_flights(Bat,batdate,logger,start_buffer,end_buffer,cluster)

    %% Plot the echolocation time stamps on the flight traces

    % load in the ciholas data + ephys data
    exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/');
    load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
    load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure_',num2str(Bat),'.mat'));
    load(strcat(exp_data_path,'ciholas/aligned_bat_position_data_',num2str(Bat),'.mat'));

    % load in confirmed echolocation events
    load(strcat(exp_data_path,'audio/saved_echolocation_timestamps.mat'));
    all_audio_events = horzcat(cell2mat(mic_audio_events));

    % Pick out echolocation event samples that are too close to 3ach other 
    all_audio_events = bootleg_prune_echo_samples(all_audio_events);
           
    flights = []; valid_flight_samples = [];
    for i=1:length(ciholas_flight_struct_resort)
        if ciholas_flight_struct_resort{i}.fclus == cluster
            flights = [flights,i];
            valid_flight_samples = [valid_flight_samples,ciholas_flight_struct_resort{i}.fstart_trimmed-start_buffer:ciholas_flight_struct_resort{i}.fend_trimmed+end_buffer];
        end
    end

    useable_echo_timestamps = [];
    for j=1:length(all_audio_events)
        if ismember(round(all_audio_events(j)*120),valid_flight_samples)
            useable_echo_timestamps = [useable_echo_timestamps,round(all_audio_events(j)*120)];
        end
    end

    % Plot in 2D.
    figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    title(strcat("Bat ",num2str(Bat)," - date - ",num2str(batdate)," - Cluster ",num2str(cluster)," flights. Echolocations = r*"));
    for i=1:length(flights)
        plot3(ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,1),...
              ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,2),...
              ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,3),'Color',[0.8 0.8 0.8]);
        scatter3(ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed,1),...
              ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed,2),...
              ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed,3),'g');
        scatter3(ciholas_r(useable_echo_timestamps,1),ciholas_r(useable_echo_timestamps,2),ciholas_r(useable_echo_timestamps,3),'*r');
    end
    hold off;
        
     % Plot in 1D. linearized
    all_flights_idxs = [];
    F_temp.flight = {}; F_temp.smp1 = []; F_temp.smp2 = []; F_temp.fclus = []; F_temp.length_m = []; F_temp.pos_stacked = [];
    F_temp.arclen = {}; F_temp.arclen_linearized = {}; F_temp.arclen_linearized_dimswap = {}; F_temp.lin_tr_block_counts = []; 
    F_temp.arclen_stacked = []; F_temp.arclen_stacked_odd = []; F_temp.arclen_stacked_even = [];
    F_temp.arclen_celled = {}; F_temp.pos_celled = {}; F_temp.arclen_celled_even = {}; F_temp.pos_celled_even = {};F_temp.arclen_celled_odd = {}; F_temp.pos_celled_odd = {};
    for jj=1:length(ciholas_flight_struct_resort)
        if ~isempty(ciholas_flight_struct_resort{jj})
            if ciholas_flight_struct_resort{jj}.fclus == cluster
                all_flights_idxs = [all_flights_idxs,jj];
                flight = ciholas_r(ciholas_flight_struct_resort{jj}.fstart_trimmed-start_buffer:ciholas_flight_struct_resort{jj}.fend_trimmed+end_buffer,:);
                audio_events_during_flight = all_audio_events
                F_temp.flight{end+1} = flight;
                F_temp.smp1 = [F_temp.smp1,ciholas_flight_struct_resort{jj}.fstart_trimmed];
                F_temp.smp2 = [F_temp.smp2,ciholas_flight_struct_resort{jj}.fend_trimmed];
                F_temp.fclus = [F_temp.fclus,ciholas_flight_struct_resort{jj}.fclus];
                flight_arclen = [];
                for k=2:length(flight)
                    px = flight(1:k,1); py = flight(1:k,2); pz = flight(1:k,3);
                    [arclen,seglen] = arclength(px,py,pz);
                    flight_arclen(k) = arclen;
                end
                F_temp.arclen{end+1} = flight_arclen;
                F_temp.length_m = [F_temp.length_m,flight_arclen(end)];
                F_temp.arclen_linearized{end+1} = flight_arclen./max(flight_arclen);
                F_temp.arclen_linearized_dimswap{end+1} = (flight_arclen./max(flight_arclen))';
                F_temp.pos_stacked = [F_temp.pos_stacked;flight(:,1:2)];
                F_temp.pos_celled{end+1} = flight(:,1:2);
            end
        end
    end
    for jj=1:length(F_temp.arclen_linearized)
        if mod(jj,2) == 1
            F_temp.arclen_stacked_odd = [F_temp.arclen_stacked_odd;F_temp.arclen_linearized{jj}'];
            F_temp.arclen_celled_odd{end+1} = F_temp.arclen_linearized{jj};
        elseif mod(jj,2) == 0
            F_temp.arclen_stacked_even = [F_temp.arclen_stacked_even;F_temp.arclen_linearized{jj}'];
            F_temp.arclen_celled_even{end+1} = F_temp.arclen_linearized{jj};
        end
        F_temp.arclen_stacked = [F_temp.arclen_stacked;F_temp.arclen_linearized{jj}'];
        F_temp.arclen_celled{end+1} = F_temp.arclen_linearized{jj}; 

        if mod(jj,2) == 1
            F_temp.pos_celled_odd{end+1} = F_temp.pos_celled{jj};
        elseif mod(jj,2) == 0
            F_temp.pos_celled_even{end+1} = F_temp.pos_celled{jj};
        end
    end

    figure(); hold on;
    scatter3(F_temp.arclen_celled{i});

end