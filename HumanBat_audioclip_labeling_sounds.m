%% Run after HumanBat_visualize_audio_event_aligned_to_flight.m

motu_Fs = 192000;

% Get best microphone for this segment 
mic_tally = zeros(1,4);
for i=1:length(above_zero_events)
    for j=1:4
        if ismember(above_zero_events(i),event_locs{j}-event_vec_start)
            mic_tally(j) = mic_tally(j)+1;
        end
    end
end

mic_use = find(mic_tally==max(mic_tally));
audio_file = dir(strcat(audio_base,'\audioConCat_mic_',num2str(mic_use(1)),'_segment_chunk_',num2str(chunk),'_break_*.mat'));
load(strcat(audio_base,audio_file.name));

% Apply segment bounds to desired segment to plot 
if seg_end>len_chunk
    signal = audioConCat(seg_start:len_chunk);
else
    signal = audioConCat(seg_start:seg_end);
end

max_seg_num = floor(length(signal)/(5*192000));
if max_seg_num==0
    max_seg_num=1;
end

for i=1:max_seg_num
    if i==max_seg_num
        seg_to_plot = signal((i-1)*5*motu_Fs+1:end);
        seg_to_plot_white = sum_whitenedSignal_allmics((i-1)*5*motu_Fs+1:end);
        aze = above_zero_events(above_zero_events<length(sum_whitenedSignal_allmics) & above_zero_events>(i-1)*5*motu_Fs+1);
        aze_plot = aze - (i-1)*5*motu_Fs+1;
    else
        seg_to_plot = signal((i-1)*5*motu_Fs+1:(i)*5*192000);
        seg_to_plot_white = sum_whitenedSignal_allmics((i-1)*5*motu_Fs+1:(i)*5*192000);
        aze = above_zero_events(above_zero_events<(i)*5*192000 & above_zero_events>(i-1)*5*motu_Fs+1);
        aze_plot = aze - (i-1)*5*motu_Fs+1;
    end

    % Plot the spectrogram with dots 
    overlap=2000;
    tscale=2;
    N=2048;
    nfft=2^nextpow2(N);
    low=2.9;
    high=10;
    
    t_=-N/2+1:N/2;
    sigma=(tscale/1e3)*motu_Fs;
    w = exp(-(t_/sigma).^2);
    dw = -2*w.*(t_/(sigma^2));
    
    % Calculate spectrogram
    [S,F,T,P]=spectrogram(seg_to_plot_white,w,overlap,nfft,motu_Fs);
    [S2]=spectrogram(seg_to_plot_white,dw,overlap,nfft,motu_Fs);
    IMAGE=100*((abs(S)+abs(S2))/2);
    IMAGE=log(IMAGE+2);
    IMAGE(IMAGE>high)=high;
    IMAGE(IMAGE<low)=low;
    IMAGE=IMAGE-low;
    IMAGE=IMAGE/(high-low);
    IMAGE=63*(IMAGE); % 63 is the cutoff for GIF
    IMAGE_MODS = log(abs(IMAGE)+1e+2);
    % Above from fb_pretty_sonogram code in FinchScope repo
    IMAGE_MODS_SUM = sum(IMAGE_MODS);

    spectrogram_ee = NaN(length(T),1);
    spectrogram_ee_idxs = floor(aze_plot/(length(seg_to_plot_white)/length(T)));
    if ~isempty(aze_plot)
        if spectrogram_ee_idxs(1) == 0; spectrogram_ee_idxs(1) = 1; end;
    end
    spectrogram_ee(spectrogram_ee_idxs) = 60000;

    % Play again?
    str = input('Play sound: (y = continue, q = quit): ', 's');
    while strcmp(str,'y')
        sound(seg_to_plot*10,motu_Fs);
        str = input('Re-play sound? (y = yes, n = no): ', 's');
        clear sound;
    end

    % Label the x_coords as audio events
    str = input('Label sounds: (c = continue, q = quit): ', 's');
    while strcmp(str,'c')
        figure(); hold on; 
        colormap(hot)
        imagesc(T,F,IMAGE_MODS); set(gca,'YDir','normal');
        axis tight;
        plot(T, spectrogram_ee','y*','MarkerSize',5); % Mark the peaks
        ylabel('Frequency kHz','FontWeight','bold');
        xlabel('time (s)','FontWeight','bold');
        title('Click audio events you think are real. Press Enter when done.','FontWeight','bold','Color','w');
        set(gca, 'Color', 'k');
        set(gca, 'XColor', 'w');
        set(gca, 'YColor', 'w');
        set(gcf, 'Color', 'k'); 
        hold off;
        [x_coords, y_coords] = ginput_bootleg;
        str = input('Re-label sounds? (y = yes, n = no): ', 's');
    end

    %% Save 50ms clips
    to_discard = [];
    for cc = 1:length(x_coords)

        % Find the aze closest to that coords 
        closest_vec = above_zero_events - (x_coords(cc)*192000)+(i-1)*5*motu_Fs+1;
        cvv = find(abs(closest_vec) == min(abs(closest_vec)));

        for jj = 1:4
            if ismember(above_zero_events(cvv),event_locs{jj} - event_vec_start)
                mic_use = j;
                audio_file = dir(strcat(audio_base,'/audioConCat_mic_',num2str(mic_use),'_segment_chunk_',num2str(chunk),'_break_*.mat'));
                load(strcat(audio_base,audio_file.name));    
                if seg_end>len_chunk
                    signal = audioConCat(seg_start:len_chunk);
                else
                    signal = audioConCat(seg_start:seg_end);
                end

                if i==5
                    seg_to_plot = signal((i-1)*5*motu_Fs+1:end);
                    seg_to_plot_white = sum_whitenedSignal_allmics((i-1)*5*motu_Fs+1:end);
                else
                    seg_to_plot = signal((i-1)*5*motu_Fs+1:(i)*5*192000);
                    seg_to_plot_white = sum_whitenedSignal_allmics((i-1)*5*motu_Fs+1:(i)*5*192000);
                end
            end
        end

        x_coords_real_sample = (x_coords(cc)*192000)+(i-1)*5*motu_Fs+1 + seg_start;
        x_coord_sample = x_coords(cc)*192000;

        % Calculate spectrogram 
        buff = 0.25;
        if x_coord_sample < buff*motu_Fs
            stp_w = seg_to_plot_white(1:round(x_coord_sample+(buff*motu_Fs)));
            stp = seg_to_plot(1:round(x_coord_sample+(buff*motu_Fs)));
        elseif x_coords(cc)*motu_Fs+(buff*motu_Fs) > length(seg_to_plot)
            stp_w = seg_to_plot_white(round(x_coord_sample-(buff*motu_Fs)):end);
            stp = seg_to_plot(round(x_coord_sample-(buff*motu_Fs)):end);
        else
            stp_w = seg_to_plot_white(round(x_coord_sample-(buff*motu_Fs)):round(x_coord_sample+(buff*motu_Fs)));
            stp = seg_to_plot(round(x_coord_sample-(buff*motu_Fs)):round(x_coord_sample+(buff*motu_Fs)));
        end

        sound(stp*10,motu_Fs);

        [S,F,T,P]=spectrogram(stp_w,w,overlap,nfft,motu_Fs);
        [S2]=spectrogram(stp_w,dw,overlap,nfft,motu_Fs);
        IMAGE=100*((abs(S)+abs(S2))/2);
        IMAGE=log(IMAGE+2);
        IMAGE(IMAGE>high)=high;
        IMAGE(IMAGE<low)=low;
        IMAGE=IMAGE-low;
        IMAGE=IMAGE/(high-low);
        IMAGE=63*(IMAGE); 
        IMAGE_MODS = log(abs(IMAGE)+1e+2);
        % Above from fb_pretty_sonogram code in FinchScope repo
        IMAGE_MODS_SUM = sum(IMAGE_MODS);
    
        figure(); hold on; 
        colormap(hot)
        imagesc(T,F,IMAGE_MODS); set(gca,'YDir','normal');
        axis tight;
        ylabel('Frequency Hz','FontWeight','bold');
        xlabel('time (s)','FontWeight','bold');
        title('50ms clip. Label!','FontWeight','bold','Color','w');
        set(gca, 'Color', 'k');
        set(gca, 'XColor', 'w');
        set(gca, 'YColor', 'w');
        set(gcf, 'Color', 'k'); 
        hold off;

        unlabeled = 1;
        unlabeled_list = {'r','rs','rc','rsrc','rcrs'};
        labeled_list = {'e','w','f','h','m','k','n'};
        str = input('Replay sound: "r". Replay slowly: "rs". Replay with more context: "rc". Replay slowly with more context: "rcrs". Label this sound: ("e" = echolocation, "w" = wingbeats, "f" = feeder, "h" = human click, "m" = madeleine voice, "k" = kevin voice, "n" = none): ', 's');
        if ismember(str,unlabeled_list)
            unlabeled = 1;
        elseif ismember(str,labeled_list)
            unlabeled = 0;
        else
            disp("Not a valid command");
        end
        
        while unlabeled 
            if strcmp(str,'r')
                sound(stp*10,192000);
            elseif strcmp(str,'rs')
                sound(stp*10,192000/2);
            elseif strcmp(str,'rc')
                buff_ = 0.5;
                if x_coord_sample < buff_*motu_Fs
                    stp_w = seg_to_plot_white(1:round(x_coord_sample+(buff_*motu_Fs)));
                    stp = seg_to_plot(1:round(x_coord_sample+(buff_*motu_Fs)));
                elseif x_coords(cc)*motu_Fs+(buff_*motu_Fs) > length(seg_to_plot)
                    stp_w = seg_to_plot_white(round(x_coord_sample-(buff_*motu_Fs)):end);
                    stp = seg_to_plot(round(x_coord_sample-(buff_*motu_Fs)):end);
                else
                    stp_w = seg_to_plot_white(round(x_coord_sample-(buff_*motu_Fs)):round(x_coord_sample+(buff*motu_Fs)));
                    stp = seg_to_plot(round(x_coord_sample-(buff_*motu_Fs)):round(x_coord_sample+(buff*motu_Fs)));
                end
                sound(stp*10,192000);
            elseif strcmp(str,'rcrs') | strcmp(str,'rsrc')
                buff_ = 0.5;
                if x_coord_sample < buff_*motu_Fs
                    stp_w = seg_to_plot_white(1:round(x_coord_sample+(buff_*motu_Fs)));
                    stp = seg_to_plot(1:round(x_coord_sample+(buff_*motu_Fs)));
                elseif x_coords(cc)*motu_Fs+(buff_*motu_Fs) > length(seg_to_plot)
                    stp_w = seg_to_plot_white(round(x_coord_sample-(buff_*motu_Fs)):end);
                    stp = seg_to_plot(round(x_coord_sample-(buff_*motu_Fs)):end);
                else
                    stp_w = seg_to_plot_white(round(x_coord_sample-(buff_*motu_Fs)):round(x_coord_sample+(buff*motu_Fs)));
                    stp = seg_to_plot(round(x_coord_sample-(buff_*motu_Fs)):round(x_coord_sample+(buff*motu_Fs)));
                end
                sound(stp*10,192000/2);
            elseif ismember(str,labeled_list)
                unlabeled = 0;
                break;
            else
                disp("Not a valid replay command or labeling command.");
            end
            % Prompt again
            str = input('Replay sound: "r". Replay slowly: "rs". Replay with more context: "rc". Replay slowly with more context: "rcrs". Label this sound: ("e" = echolocation, "w" = wingbeats, "f" = feeder, "h" = human click, "m" = madeleine voice, "k" = kevin voice, "n" = none): ', 's');   
        end
        
        % Once labeled...
        if strcmp(str,'n')
            to_discard = [to_discard,cc];
        elseif strcmp(str,'e')
            name = strcat('echolocation_',num2str(batdate),'_mic_',num2str(mic_use),'_chunk_',num2str(chunk),'_ss_',num2str(round(x_coords_real_sample)),'.wav');
        elseif strcmp(str,'f')
            name = strcat('feederclick_',num2str(batdate),'_mic_',num2str(mic_use),'_chunk_',num2str(chunk),'_ss_',num2str(round(x_coords_real_sample)),'.wav');
        elseif strcmp(str,'h')
            name = strcat('humanclick_',num2str(batdate),'_mic_',num2str(mic_use),'_chunk_',num2str(chunk),'_ss_',num2str(round(x_coords_real_sample)),'.wav');
        elseif strcmp(str,'m')
            name = strcat('madeleinevoice_',num2str(batdate),'_mic_',num2str(mic_use),'_chunk_',num2str(chunk),'_ss_',num2str(round(x_coords_real_sample)),'.wav');
        elseif strcmp(str,'k')
            name = strcat('kevinvoice_',num2str(batdate),'_mic_',num2str(mic_use),'_chunk_',num2str(chunk),'_ss_',num2str(round(x_coords_real_sample)),'.wav');
        elseif strcmp(str,'w')
            name = strcat('wingbeats_',num2str(batdate),'_mic_',num2str(mic_use),'_chunk_',num2str(chunk),'_ss_',num2str(round(x_coords_real_sample)),'.wav');
        end

        if strcmp(str,'n')
            continue
        else
            audiowrite(strcat('C:\Users\YartsevLabComputer5\Box\Audio_Project_2024\221126\snippits_of_audio\hand_selected_snippits\',name),stp,motu_Fs);
        end

        close all;
    end
end

disp("Finished with segment! Rerun HumanBat_visualize_audio_event_aligned_to_flight.m to select a new segment to label");