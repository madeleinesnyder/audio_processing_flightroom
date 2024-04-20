function [] = bootleg_plot_echos_MvsK(Bat,batdate,logger,unit,start_buffer,end_buffer,cluster,plot_units)

    %% Plot the echolocation time stamps on the flight traces

    % load in the ciholas data + ephys data
    exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/');
    load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));
    load(strcat(exp_data_path,'ciholas/pruned_resorted_trimmed_ciholas_bat_final_flight_structure_',num2str(Bat),'.mat'));
    load(strcat(exp_data_path,'ciholas/aligned_bat_position_data_',num2str(Bat),'.mat'));
    obat = load(strcat(exp_data_path,'ciholas/aligned_bat_position_data_',num2str(Bat),'.mat'));
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

    % load in confirmed echolocation events
    load(strcat(exp_data_path,'audio/saved_echolocation_timestamps.mat'));
    all_audio_events = horzcat(cell2mat(mic_audio_events));

    % Pick out echolocation event samples that are too close to 3ach other 
    all_audio_events = bootleg_prune_echo_samples(all_audio_events);
           
    flights_K = []; valid_flight_samples_K = []; flights_M = []; valid_flight_samples_M = [];
    for i=1:length(ciholas_flight_struct_resort)
        if ciholas_flight_struct_resort{i}.fclus == cluster & ismember(i,Flight_Group_Matrix_To_K)
            flights_K = [flights_K,i];
            valid_flight_samples_K = [valid_flight_samples_K,ciholas_flight_struct_resort{i}.fstart_trimmed-start_buffer:ciholas_flight_struct_resort{i}.fend_trimmed+end_buffer];
        elseif ciholas_flight_struct_resort{i}.fclus == cluster & ismember(i,Flight_Group_Matrix_To_M)
            flights_M = [flights_M,i];
            valid_flight_samples_M = [valid_flight_samples_M,ciholas_flight_struct_resort{i}.fstart_trimmed-start_buffer:ciholas_flight_struct_resort{i}.fend_trimmed+end_buffer];
        end
    end

    valid_flight_samples = [valid_flight_samples_M,valid_flight_samples_K];

    [useable_echo_timestamps,who_echod,probably_not_echos,which_mics] = bootleg_who_echolocated(ciholas_r,obat.ciholas_r,mic_audio_events,valid_flight_samples);
    useable_echo_timestamps(probably_not_echos) = [];
    who_echod(probably_not_echos) = [];
    which_mics(probably_not_echos) = [];
    useable_echo_timestamps_bat1 = useable_echo_timestamps(who_echod==1);
    useable_echo_timestamps_bat2 = useable_echo_timestamps(who_echod==2);

    useable_echo_timestamps_M = []; useable_echo_timestamps_K = [];
    for j=1:length(useable_echo_timestamps_bat1)
        if ismember(useable_echo_timestamps_bat1(j),valid_flight_samples_K)
            useable_echo_timestamps_K = [useable_echo_timestamps_K,useable_echo_timestamps_bat1(i)];
        elseif ismember(useable_echo_timestamps_bat1(i),valid_flight_samples_M)
            useable_echo_timestamps_M = [useable_echo_timestamps_M,useable_echo_timestamps_bat1(i)];
        end
    end

    % Plot in 2D.
    figure(); 
    subplot(1,2,1); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    title(strcat("Bat ",num2str(Bat)," - date - ",num2str(batdate)," - M. Cluster ",num2str(cluster)," flights. Echolocations = r*"));
    for i=1:length(flights_M)
        plot3(ciholas_r(ciholas_flight_struct_resort{flights_M(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights_M(i)}.fend_trimmed,1),...
              ciholas_r(ciholas_flight_struct_resort{flights_M(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights_M(i)}.fend_trimmed,2),...
              ciholas_r(ciholas_flight_struct_resort{flights_M(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights_M(i)}.fend_trimmed,3),'Color',[0.8 0.8 0.8]);
        scatter3(ciholas_r(ciholas_flight_struct_resort{flights_M(i)}.fstart_trimmed,1),...
              ciholas_r(ciholas_flight_struct_resort{flights_M(i)}.fstart_trimmed,2),...
              ciholas_r(ciholas_flight_struct_resort{flights_M(i)}.fstart_trimmed,3),'g');
        scatter3(ciholas_r(useable_echo_timestamps_M,1),ciholas_r(useable_echo_timestamps_M,2),ciholas_r(useable_echo_timestamps_M,3),'*r');
    end
    hold off;
    subplot(1,2,2); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    title(strcat("Bat ",num2str(Bat)," - date - ",num2str(batdate)," - K. Cluster ",num2str(cluster)," flights. Echolocations = r*"));
    for i=1:length(flights_K)
        plot3(ciholas_r(ciholas_flight_struct_resort{flights_K(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights_K(i)}.fend_trimmed,1),...
              ciholas_r(ciholas_flight_struct_resort{flights_K(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights_K(i)}.fend_trimmed,2),...
              ciholas_r(ciholas_flight_struct_resort{flights_K(i)}.fstart_trimmed:ciholas_flight_struct_resort{flights_K(i)}.fend_trimmed,3),'Color',[0.8 0.8 0.8]);
        scatter3(ciholas_r(ciholas_flight_struct_resort{flights_K(i)}.fstart_trimmed,1),...
              ciholas_r(ciholas_flight_struct_resort{flights_K(i)}.fstart_trimmed,2),...
              ciholas_r(ciholas_flight_struct_resort{flights_K(i)}.fstart_trimmed,3),'g');
        scatter3(ciholas_r(useable_echo_timestamps_K,1),ciholas_r(useable_echo_timestamps_K,2),ciholas_r(useable_echo_timestamps_K,3),'*r');
    end
    hold off;    
        
    % Plot in 1D linearized with ephys.
    all_audio_events = all_audio_events*120;
    all_flights_idxs = [];
    F_temp.flight = {}; F_temp.smp1 = []; F_temp.smp2 = []; F_temp.fclus = []; F_temp.length_m = []; F_temp.pos_stacked = [];
    F_temp.arclen = {}; F_temp.arclen_linearized = {}; F_temp.arclen_linearized_dimswap = {}; F_temp.lin_tr_block_counts = []; 
    F_temp.arclen_stacked = []; F_temp.arclen_stacked_odd = []; F_temp.arclen_stacked_even = [];
    F_temp.arclen_celled = {}; F_temp.pos_celled = {}; F_temp.arclen_celled_even = {}; F_temp.pos_celled_even = {};F_temp.arclen_celled_odd = {}; F_temp.pos_celled_odd = {};
    F_temp.audio_during_flight = {}; F_temp.ephys_during_flight = {};
    human_ = {};
    for jj=1:length(ciholas_flight_struct_resort)
        if ~isempty(ciholas_flight_struct_resort{jj})
            if ciholas_flight_struct_resort{jj}.fclus == cluster
                nearest_ephys_points = [];
                all_flights_idxs = [all_flights_idxs,jj];
                flight = ciholas_r(ciholas_flight_struct_resort{jj}.fstart_trimmed-start_buffer:ciholas_flight_struct_resort{jj}.fend_trimmed+end_buffer,:);
                audio_events_during_flight = all_audio_events(all_audio_events >= (ciholas_flight_struct_resort{jj}.fstart_trimmed-start_buffer) & all_audio_events <= (ciholas_flight_struct_resort{jj}.fend_trimmed+end_buffer));
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
                if ~isempty(audio_events_during_flight)
                    F_temp.audio_during_flight{end+1} = audio_events_during_flight-(ciholas_flight_struct_resort{jj}.fstart_trimmed-start_buffer);
                else
                    F_temp.audio_during_flight{end+1} = audio_events_during_flight;
                end
                if ismember(jj,Flight_Group_Matrix_To_M)
                    human_{end+1} = 1;
                elseif ismember(jj,Flight_Group_Matrix_To_K)
                    human_{end+1} = 2;
                else
                    human_{end+1} = 0;
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

    % Bin and convolve signal with gaussian kernel
    bins = [0:0.001:1]; bin_counts_K = zeros(1,length(bins)-1); bin_counts_M = zeros(1,length(bins)-1);
    for jj=1:length(F_temp.audio_during_flight)
        if human_{jj} == 2
            full_vec = [1,F_temp.audio_during_flight{jj},length(F_temp.flight{jj})];
            scaled_audio_times = (full_vec - min(full_vec)) / (max(full_vec) - min(full_vec));
            scaled_audio_times([1,end]) = [];
            for ss = 1:length(scaled_audio_times)
                bin_dist = bins-scaled_audio_times(ss);
                pos_bin_dist = bin_dist(bin_dist>0);
                [minval,minpos] = min(pos_bin_dist);
                bin_add = find(bin_dist == minval, 1, 'first');
                bin_counts_K(bin_add) = bin_counts_K(bin_add)+1;
            end
        elseif human_{jj} == 1
            full_vec = [1,F_temp.audio_during_flight{jj},length(F_temp.flight{jj})];
            scaled_audio_times = (full_vec - min(full_vec)) / (max(full_vec) - min(full_vec));
            scaled_audio_times([1,end]) = [];
            for ss = 1:length(scaled_audio_times)
                bin_dist = bins-scaled_audio_times(ss);
                pos_bin_dist = bin_dist(bin_dist>0);
                [minval,minpos] = min(pos_bin_dist);
                bin_add = find(bin_dist == minval, 1, 'first');
                bin_counts_M(bin_add) = bin_counts_M(bin_add)+1;
            end
        end
    end

    % Define gaussian kernel
    sigma = 1;
    N = length(bin_counts_M)/10-1; 
    x = linspace(-3*sigma,3*sigma, N);
    gaussianKernel = exp(-x.^2 / (2*sigma^2));
    gaussianKernel = gaussianKernel / sum(gaussianKernel);
    
    bin_counts_conv_M = conv(bin_counts_M,gaussianKernel,'same');
    bin_counts_conv_K = conv(bin_counts_K,gaussianKernel,'same');

    figure(); 
    subplot(1,2,1); hold on; xlim([0 1000]); title(strcat("All cluster ",num2str(cluster)," flights to Madeleine"));  
    ctr = 0;
    for jj=1:length(F_temp.audio_during_flight)
        if human_{jj} == 1
            ctr = ctr+0.1;
            plot([0 1000], [ctr ctr],'Color',[0.8 0.8 0.8],'LineWidth',0.2); 
            if ~isempty(F_temp.audio_during_flight{jj})
                full_vec = [1,F_temp.audio_during_flight{jj},length(F_temp.flight{jj})];
                % Scale between 0 and 1 and plot
                scaled_audio_times = (full_vec - min(full_vec)) / (max(full_vec) - min(full_vec));
                scaled_audio_times([1,end]) = [];
                y_values = 0.5* ones(size(scaled_audio_times));
                plot(scaled_audio_times*1000,ctr,'r*','MarkerSize',5,'LineWidth',2);
            end
        end
    end
    plot(bin_counts_conv_M,'g');
    hold off;

    subplot(1,2,2); hold on; xlim([0 1000]); title(strcat("All cluster ",num2str(cluster)," flights to Kevin")); 
    ctr = 0;
    for jj=1:length(F_temp.audio_during_flight)
        if human_{jj} == 2
            ctr = ctr+0.1;
            plot([0 1000], [ctr ctr],'Color',[0.8 0.8 0.8],'LineWidth',0.2); 
            if ~isempty(F_temp.audio_during_flight{jj})
                full_vec = [1,F_temp.audio_during_flight{jj},length(F_temp.flight{jj})];
                % Scale between 0 and 1 and plot
                scaled_audio_times = (full_vec - min(full_vec)) / (max(full_vec) - min(full_vec));
                scaled_audio_times([1,end]) = [];
                y_values = 0.5* ones(size(scaled_audio_times));
                plot(scaled_audio_times*1000,ctr,'r*','MarkerSize',5,'LineWidth',2);
            end
        end
    end
    plot(bin_counts_conv_K,'r');
    hold off;

    %% Bayesian approach modeling an echoloation train
    % Define prior (Beta distribution)
%     alpha_prior = ones(1,1000);
%     beta_prior = ones(1,1000);
% 
%     alpha_post = alpha_prior + bin_counts_M;
%     beta_post = beta_prior + repmat(length(flights_M),1,1000) - bin_counts_M;
% 
%     post_mean = alpha_post ./ (alpha_post + beta_post);
%     post_var = (alpha_post .* beta_post) ./ ((alpha_post + beta_post).^2 .* (alpha_post + beta_post + 1));
% 
%     new_samples = betarnd(mean(alpha_post),mean(beta_post),[1000,1]);
    
end