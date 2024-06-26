function [M_flight_hz_samples,K_flight_hz_samples,M_flight_hz_samples_T,K_flight_hz_samples_T] = bootleg_plot_echos_on_flights(Bat,batdate,logger,unit,start_buffer,end_buffer,cluster,plot_units)

    %% Plot the echolocation time stamps on the flight traces

    % load in the ciholas data + ephys data
    exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/');
    load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
    load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure_',num2str(Bat),'.mat'));
    load(strcat(exp_data_path,'ciholas/aligned_bat_position_data_',num2str(Bat),'.mat'));

    obat = load(strcat(exp_data_path,'ciholas/aligned_bat_position_data_14643.mat'));

    % load in CONFIRMED echolocation events (these are NOT confirmed
    % currently)
    load(strcat(exp_data_path,'audio/saved_echolocation_timestamps.mat'));
    all_audio_events = horzcat(cell2mat(mic_audio_events));

    % Pick out echolocation event samples that are too close to 3ach other 
    all_audio_events = bootleg_prune_echo_samples(all_audio_events);
           
    % For all flights in a cluster, extract the samples in all these
    % flights
    flights = []; valid_flight_samples = [];
    for i=1:length(ciholas_flight_struct_resort)
        if ciholas_flight_struct_resort{i}.fclus == cluster
            flights = [flights,i];
            valid_flight_samples = [valid_flight_samples,ciholas_flight_struct_resort{i}.fstart_trimmed-start_buffer:ciholas_flight_struct_resort{i}.fend_trimmed+end_buffer];
        end
    end

    [useable_echo_timestamps,who_echod,probably_not_echos,which_mics] = bootleg_who_echolocated(ciholas_r,obat.ciholas_r,mic_audio_events,valid_flight_samples);
    useable_echo_timestamps(probably_not_echos) = [];
    who_echod(probably_not_echos) = [];
    which_mics(probably_not_echos) = [];
    useable_echo_timestamps_bat1 = useable_echo_timestamps(who_echod==1);
    useable_echo_timestamps_bat2 = useable_echo_timestamps(who_echod==2);

    %% Plot in 2D.
    figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    title(strcat("Bat ",num2str(Bat)," - date - ",num2str(batdate)," - Cluster ",num2str(cluster)," flights. Echolocations = r*"));
    for i=1:length(flights)
        plot3(ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,1),...
              ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,2),...
              ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,3),'Color',[0.8 0.8 0.8]);
        plot3(obat.ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,1),...
              obat.ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,2),...
              obat.ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,3),'Color',[0.8 1 0.8]);
    end
    for i=1:length(useable_echo_timestamps)
        if who_echod(i) == 1
            scatter3(ciholas_r(useable_echo_timestamps(i),1),ciholas_r(useable_echo_timestamps(i),2),ciholas_r(useable_echo_timestamps(i),3),'*r');
        elseif who_echod(i) == 2
            scatter3(obat.ciholas_r(useable_echo_timestamps(i),1),obat.ciholas_r(useable_echo_timestamps(i),2),obat.ciholas_r(useable_echo_timestamps(i),3),'*k');
        end
    end
    hold off;

    % Plot in 2D ONLY BAT 1
    figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    title(strcat("Bat ",num2str(Bat)," - date - ",num2str(batdate)," - Cluster ",num2str(cluster)," flights. Echolocations = r*"));
    for i=1:length(flights)
        plot3(ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,1),...
              ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,2),...
              ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,3),'Color',[0.8 0.8 0.8]);
    end
    for i=1:length(useable_echo_timestamps)
        if who_echod(i) == 1
            scatter3(ciholas_r(useable_echo_timestamps(i),1),ciholas_r(useable_echo_timestamps(i),2),ciholas_r(useable_echo_timestamps(i),3),'*r');
        end
    end
    hold off;
        

    %% For humans split
    load(strcat(exp_data_path,'ciholas/aligned_human_position_data','.mat'));
    madeleine=2; kevin=4;

    Flight_Group_Matrix_To_M=[]; Flight_Group_Matrix_To_K=[];  Flight_Group_Matrix_From_K=[];  Flight_Group_Matrix_From_M=[]; 
    for i=1:length(ciholas_flight_struct_resort)
        if ~isempty(ciholas_flight_struct_resort{i}) & ciholas_flight_struct_resort{i}.fclus == cluster
            % Find if flight was to Madeleine
            mean_end_position = mean(ciholas_r(ciholas_flight_struct_resort{i}.fend_trimmed-50:ciholas_flight_struct_resort{i}.fend_trimmed,:));
            mean_madeleine_position = mean(human_r(ciholas_flight_struct_resort{i}.fend_trimmed-50:ciholas_flight_struct_resort{i}.fend_trimmed,:,madeleine));
            mean_kevin_position = mean(human_r(ciholas_flight_struct_resort{i}.fend_trimmed-50:ciholas_flight_struct_resort{i}.fend_trimmed,:,kevin));
            if pdist2(mean_end_position, mean_madeleine_position, 'euclidean') < 600
                Flight_Group_Matrix_To_M = [Flight_Group_Matrix_To_M,i];
            elseif pdist2(mean_end_position, mean_kevin_position, 'euclidean') < 600
                Flight_Group_Matrix_To_K = [Flight_Group_Matrix_To_K,i];
            end
            mean_start_position =mean(ciholas_r(ciholas_flight_struct_resort{i}.fstart_trimmed-50:ciholas_flight_struct_resort{i}.fstart_trimmed,:));
            mean_madeleine_position = mean(human_r(ciholas_flight_struct_resort{i}.fstart_trimmed-50:ciholas_flight_struct_resort{i}.fstart_trimmed,:,madeleine));
            mean_kevin_position = mean(human_r(ciholas_flight_struct_resort{i}.fstart_trimmed-50:ciholas_flight_struct_resort{i}.fstart_trimmed,:,kevin));
            if pdist2(mean_start_position, mean_madeleine_position, 'euclidean') < 600
                Flight_Group_Matrix_From_M = [Flight_Group_Matrix_From_M,i];
            elseif pdist2(mean_start_position, mean_kevin_position, 'euclidean') < 600
                Flight_Group_Matrix_From_K = [Flight_Group_Matrix_From_K,i];
            end
        end
    end

    % Plot in 2D ONLY BAT 1
    figure(); 
    subplot(1,2,1); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    title(strcat("Bat ",num2str(Bat)," - date - ",num2str(batdate)," - Cluster ",num2str(cluster)," flights to M. Echolocations = r*"));
    M_samples = []; M_hz_window_samples = {}; M_hz_window_samples_T = {};
    for i=1:length(flights)
        if ismember(flights(i),Flight_Group_Matrix_To_M)
            M_samples = [M_samples,ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed];
            M_hz_window_samples{end+1} = [ciholas_flight_struct_resort{flights(i)}.fend_trimmed-120:ciholas_flight_struct_resort{flights(i)}.fend_trimmed];
            M_hz_window_samples_T{end+1} = [ciholas_flight_struct_resort{flights(i)}.fstart_trimmed-60:ciholas_flight_struct_resort{flights(i)}.fstart_trimmed+60];
            plot3(ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,1),...
                  ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,2),...
                  ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,3),'Color',[0.8 0.8 0.8]);
        
        end
    end
    for i=1:length(useable_echo_timestamps)
        if who_echod(i) == 1 & ismember(useable_echo_timestamps(i),M_samples)
            scatter3(ciholas_r(useable_echo_timestamps(i),1),ciholas_r(useable_echo_timestamps(i),2),ciholas_r(useable_echo_timestamps(i),3),'*r');
        end
    end
    hold off;
    subplot(1,2,2); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    title(strcat("Bat ",num2str(Bat)," - date - ",num2str(batdate)," - Cluster ",num2str(cluster)," flights to K. Echolocations = r*"));
    K_samples = []; K_hz_window_samples = {}; K_hz_window_samples_T = {};
    for i=1:length(flights)
        if ismember(flights(i),Flight_Group_Matrix_To_K)
            K_samples = [K_samples,ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed];
            K_hz_window_samples{end+1} = [ciholas_flight_struct_resort{flights(i)}.fend_trimmed-120:ciholas_flight_struct_resort{flights(i)}.fend_trimmed];
            K_hz_window_samples_T{end+1} = [ciholas_flight_struct_resort{flights(i)}.fstart_trimmed-60:ciholas_flight_struct_resort{flights(i)}.fstart_trimmed+60];
            plot3(ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,1),...
                  ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,2),...
                  ciholas_r(ciholas_flight_struct_resort{flights(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights(i)}.fend_trimmed,3),'Color',[0.8 0.8 0.8]);
        end
    end
    for i=1:length(useable_echo_timestamps)
        if who_echod(i) == 1 & ismember(useable_echo_timestamps(i),K_samples)
            scatter3(ciholas_r(useable_echo_timestamps(i),1),ciholas_r(useable_echo_timestamps(i),2),ciholas_r(useable_echo_timestamps(i),3),'*r');
        end
    end
    hold off;

    % Get echolocation rate in last 1 second before landing
    M_flight_hz_samples = [];
    for i=1:length(M_hz_window_samples)
        if isempty(M_hz_window_samples{i})
            M_flight_hz_samples = [M_flight_hz_samples,0];
        else
            temp_useable_echos = [];
            for j=1:length(useable_echo_timestamps)
                if ismember(useable_echo_timestamps(j),M_hz_window_samples{i})
                    temp_useable_echos = [temp_useable_echos,useable_echo_timestamps(j)];
                end
            end
            M_flight_hz_samples = [M_flight_hz_samples,length(temp_useable_echos)];
        end
    end

    K_flight_hz_samples = [];
    for i=1:length(K_hz_window_samples)
        if isempty(K_hz_window_samples{i})
            K_flight_hz_samples = [K_flight_hz_samples,0];
        else
            temp_useable_echos = [];
            for j=1:length(useable_echo_timestamps)
                if ismember(useable_echo_timestamps(j),K_hz_window_samples{i})
                    temp_useable_echos = [temp_useable_echos,useable_echo_timestamps(j)];
                end
            end
            K_flight_hz_samples = [K_flight_hz_samples,length(temp_useable_echos)];
        end
    end
    % Get echolocation rate in last 1 second before landing
    M_flight_hz_samples_T = [];
    for i=1:length(M_hz_window_samples_T)
        if isempty(M_hz_window_samples_T{i})
            M_flight_hz_samples_T = [M_flight_hz_samples_T,0];
        else
            temp_useable_echos = [];
            for j=1:length(useable_echo_timestamps)
                if ismember(useable_echo_timestamps(j),M_hz_window_samples_T{i})
                    temp_useable_echos = [temp_useable_echos,useable_echo_timestamps(j)];
                end
            end
            M_flight_hz_samples_T = [M_flight_hz_samples_T,length(temp_useable_echos)];
        end
    end

    K_flight_hz_samples_T = [];
    for i=1:length(K_hz_window_samples_T)
        if isempty(K_hz_window_samples_T{i})
            K_flight_hz_samples_T = [K_flight_hz_samples_T,0];
        else
            temp_useable_echos = [];
            for j=1:length(useable_echo_timestamps)
                if ismember(useable_echo_timestamps(j),K_hz_window_samples_T{i})
                    temp_useable_echos = [temp_useable_echos,useable_echo_timestamps(j)];
                end
            end
            K_flight_hz_samples_T = [K_flight_hz_samples_T,length(temp_useable_echos)];
        end
    end

    %% Plot in 1D linearized with ephys.
    all_audio_events = all_audio_events*120;
    all_flights_idxs = [];
    F_temp.flight = {}; F_temp.smp1 = []; F_temp.smp2 = []; F_temp.fclus = []; F_temp.length_m = []; F_temp.pos_stacked = [];
    F_temp.arclen = {}; F_temp.arclen_linearized = {}; F_temp.arclen_linearized_dimswap = {}; F_temp.lin_tr_block_counts = []; 
    F_temp.arclen_stacked = []; F_temp.arclen_stacked_odd = []; F_temp.arclen_stacked_even = [];
    F_temp.arclen_celled = {}; F_temp.pos_celled = {}; F_temp.arclen_celled_even = {}; F_temp.pos_celled_even = {};F_temp.arclen_celled_odd = {}; F_temp.pos_celled_odd = {};
    F_temp.audio_during_flight = {}; F_temp.ephys_during_flight = {};
    for jj=1:length(ciholas_flight_struct_resort)
        if ~isempty(ciholas_flight_struct_resort{jj})
            if ciholas_flight_struct_resort{jj}.fclus == cluster
                nearest_ephys_points = [];
                all_flights_idxs = [all_flights_idxs,jj];
                flight = ciholas_r(ciholas_flight_struct_resort{jj}.fstart_trimmed-start_buffer:ciholas_flight_struct_resort{jj}.fend_trimmed+end_buffer,:);
                audio_events_during_flight_bat1 = useable_echo_timestamps_bat1(useable_echo_timestamps_bat1 >= (ciholas_flight_struct_resort{jj}.fstart_trimmed-start_buffer) & useable_echo_timestamps_bat1 <= (ciholas_flight_struct_resort{jj}.fend_trimmed+end_buffer));
                flight_dur = (ciholas_flight_struct_resort{jj}.fend_trimmed -  ciholas_flight_struct_resort{jj}.fstart_trimmed + start_buffer + end_buffer);
                flight_vec = [1:flight_dur]./120;
                if ciholas_flight_struct_resort{jj}.ephys_trimmed{unit} == 0
                    disp(strcat("No unit firing"));
                else
                    for ul=1:length(ciholas_flight_struct_resort{jj}.ephys_trimmed{unit})
                        nearest_ephys_points = [nearest_ephys_points,dsearchn(flight_vec',ciholas_flight_struct_resort{jj}.ephys_trimmed{unit}(ul)/1e6)];
                    end
                end
                F_temp.ephys_during_flight{end+1} = nearest_ephys_points;
                F_temp.flight{end+1} = flight;
                F_temp.smp1 = [F_temp.smp1,ciholas_flight_struct_resort{jj}.fstart_trimmed-start_buffer];
                F_temp.smp2 = [F_temp.smp2,ciholas_flight_struct_resort{jj}.fend_trimmed+end_buffer];
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
                if ~isempty(audio_events_during_flight_bat1)
                    F_temp.audio_during_flight{end+1} = audio_events_during_flight_bat1-(ciholas_flight_struct_resort{jj}.fstart_trimmed-start_buffer);
                else
                    F_temp.audio_during_flight{end+1} = audio_events_during_flight_bat1;
                end
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

    figure(); hold on; xlim([0 1]);  
    ctr = 0;
    for jj=1:length(F_temp.audio_during_flight)
        ctr = ctr+0.1;
        plot([0 1], [ctr ctr],'k-','LineWidth',0.2); 
        if ~isempty(F_temp.audio_during_flight{jj})
            full_vec = [1,F_temp.audio_during_flight{jj},length(F_temp.flight{jj})];
            % Scale between 0 and 1 and plot
            scaled_audio_times = (full_vec - min(full_vec)) / (max(full_vec) - min(full_vec));
            scaled_audio_times([1,end]) = [];
            y_values = 0.5* ones(size(scaled_audio_times));
            plot(scaled_audio_times,ctr,'r*','MarkerSize',5,'LineWidth',2);
        end
    end
    hold off;

end