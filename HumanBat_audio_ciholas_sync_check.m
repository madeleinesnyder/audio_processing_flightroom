
% Not finished, weird stuff with the alignment plots 
% Run after HumanBat_get_ttl_peaks.m

% Sampling rates
motu_Fs = 192000; % audio sampling rate
ciholas_Fs = 120; % ciholas sampling rate

Bat = 32626; logger = 15;
alignment_ = 'start';

exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/');

% Load in the flight data 
load(strcat(exp_data_path,'ciholas/aligned_bat_position_data_',num2str(Bat),'.mat'));

% Load in the Extracted Behavior data
load(strcat(exp_data_path,'ciholas/Extracted_Behavior_', num2str(batdate),'_',num2str(Bat),'.mat'));

% Load in spike data
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));


% add in the last sample for 221125 because there is one TTL in the final
% file that wasn't concatenated.
%pk_loc_vec = [pk_loc_vec;pk_loc_vec(end)+(pk_loc_vec(end)-pk_loc_vec(end-1))];

if strcmp(alignment_,'ending')
    % When aligning using the LAST ttl, 
    % Add buffer to the start. The buffer is how much longer the ciholas_r vector is compared to the audio ttl vector
    gap_to_fill = (length(ciholas_r)/ciholas_Fs) - (last_ttl_sample/motu_Fs);
    % Scale this new gap by the audio sampling rate (motu samping rate)
    gap_to_fill_motu = gap_to_fill*motu_Fs;
    
    % Add this buffer to the ttl vector 
    pk_loc_vec_zerod = pk_loc_vec + gap_to_fill_motu;
    % Scale this ttl vector to the ciholas sample rate 
    pk_loc_vec_zerod_cih = round(pk_loc_vec_zerod / motu_Fs * ciholas_Fs);
    % Make a binary vector Y for all ttl times. 
    Y = NaN(1,pk_loc_vec_zerod_cih(end));
    Y(pk_loc_vec_zerod_cih) = 1;

    %% Load in the sorted_event times and add the same amount
    sorted_events_zerod = sorted_events + gap_to_fill_motu;
    % Sample the sorted_events_zerod to seconds, and then to samples 
    sorted_events_zerod_cih = round(sorted_events_zerod / motu_Fs * ciholas_Fs);
    % Make a binary vector Z for all sorted event times. 
    Z = NaN(1,pk_loc_vec_zerod_cih(end));
    Z(sorted_events_zerod_cih) = 1;

elseif strcmp(alignment_,'start')

    % If you want to align to the START... subtract the first ttl from all values.
    
    pk_loc_vec_zerod_to_start = pk_loc_vec - first_ttl_sample+1;
    pk_loc_vec_zerod_to_start_cih = round(pk_loc_vec_zerod_to_start / motu_Fs * ciholas_Fs);
    Y = NaN(1,pk_loc_vec_zerod_to_start_cih(end));
    pk_loc_vec_zerod_to_start_cih(1) = 1;
    Y(pk_loc_vec_zerod_to_start_cih) = 1;

    %% Load in the sorted_event times and add the same amount
    sorted_events_zerod_to_start = sorted_events - first_ttl_sample + 1;
    % Sample the sorted_events_zerod to seconds, and then to samples 
    sorted_events_zerod_to_start_cih = round(sorted_events_zerod_to_start / motu_Fs * ciholas_Fs);
    % Make a binary vector Z for all sorted event times. 
    Z = NaN(1,pk_loc_vec_zerod_to_start_cih(end));
    Z = NaN(1,length(ciholas_r));
    Z(sorted_events_zerod_to_start_cih) = 1;
end

% Plot the ciholas velocity, the audio event times, and the ttl times. 
figure(); 
subplot(3,1,1);
plot(t(1:end-1),diff(ciholas_r));
subplot(3,1,2);
scatter(t,Z);
subplot(3,1,3);
scatter(t,Y);

%% Now look at the wingbeats to see if the alignment was proper

% Load in the Extracted Behavior data
%load(strcat(exp_data_path,'ciholas/Extracted_Behavior_', num2str(batdate),'_',num2str(Bat),'.mat'));

% Load in spike data
%load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));

% First ephys TTL time
usec_of_first_ephys_TTL = (B_ephys_data.TT_unit(1).Timestamps(1)-B_ephys_data.TT_unit(1).AlignedTimestamps(1));
% ^^ This is the timestamp usec you want to subtract from all values in the CSC file to get the same alignment to the ephys spike file

% Load in the CSC.mat data
%csc_files = dir(fullfile(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/*CSC*')));
%load(strcat(csc_files(1).folder,'/',csc_files(1).name));
Fs = Estimated_channelFS_Transceiver(1);
raw_V = double(AD_count_int16*AD_count_to_uV_factor);
tstamps = Timestamps_of_first_samples_usec(1)+[0:length(raw_V)-1]/Fs*1e6;
tstamps_shift = tstamps-usec_of_first_ephys_TTL;   
t_stamps_g0 = find(tstamps_shift > 0); t_g0 = t_stamps_g0(1); t_stamps_g0=[];
tstamps = []; 
raw_V_shift = raw_V(t_g0:end); raw_V = [];

% % Plot the ciholas acceleration value for seg_amt samples ciholas seconds is /120
% % ephys seconds is /1e6
% seg_amt_sec = 780; % Take first 780 seconds of the session
% ciholas_sample_amt = seg_amt_sec*120;
% ephys_sample_amt = seg_amt_sec*1e6;

% Downsample the ephys data
t_sample_length = t(end)/length(t);
extra_time = tstamps_shift(end)/1e6-t(end);
extra_samples = round(extra_time/t_sample_length);
desired_ephys_trail_sample_length = length(t)+extra_samples;

oL = length(raw_V_shift);
y = interp1(1:oL,raw_V_shift,linspace(1,oL,desired_ephys_trail_sample_length)); y = y/3000;

% Plot ciholas data, wingbeats, and ephys!
Fs = 120;
figure(); hold on;
[up,lo] = envelope(a_flt+normrnd(0,1e-3,length(a_flt),1),Fs/10,'peak');     % Envelope of the acceleration signal (noise addedd to avoid problems with splining)
env = normalize(up - lo,'range');                                           % Amplitude of the envelope
env_th = otsuthresh(histcounts(env));                                       % Threshold (based on Otsu method). Can be set at 0.35
wBeats = movsum(env>env_th,1.67*Fs)>Fs/6; 
%figure(); hold on; %Heuristic criterion for flight detection
area(t,wBeats*3,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(t,normalize(v_abs,'range'));
plot(t,r(:,1),t,r(:,2));
plot(t,normalize(a_flt,'range',[-1 1]));
plot(t,normalize(movsum(env>env_th,1.67*Fs),'range'));
plot(t,up);
plot(t,lo);
plot(t,y(1:length(t)));
scatter(t,Z,'*r');
[y_up,y_lo] = envelope(y(1:length(t)),Fs/10,'peak');     % Envelope of the acceleration signal (noise addedd to avoid problems with splining)
y_env = normalize(y_up - y_lo,'range'); 
plot(t,y_up);
plot(t,y_lo);
hold off;

% prompt = "Does this plot look aligned? (Y) or (N)";
% g = input(prompt,'s');
% if g=="Y"
%     wingbeats_synced=1;
% else
%     wingbeats_synced=0;
% end

% % Calculate autocorreclaion between acceleration and ephys envelopes
% start_xcorr_rng = [find(t==1159):find(t==1414)]; end_xcorr_rng = [find(t==5744):find(t==5825)];
% [flight_ephys_xcorr,lags] = xcorr(y_up(start_xcorr_rng(1):start_xcorr_rng(end)),up(start_xcorr_rng(1):start_xcorr_rng(end)),10,'biased');
% [flight_ephys_xcorr,lags] = xcorr(y(1:length(t)),a_flt);
% figure(); stem(lags,flight_ephys_xcorr);

end





% Code to align the wingbeats of ciholas to the wingbeats in the LFP to
% double check that everything is sync'd
%exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/');

% Load in the ciholas data
load(strcat(exp_data_path,'ciholas/aligned_bat_position_data_',num2str(Bat),'.mat'));
load(strcat(exp_data_path,'ciholas/Extracted_Behavior_', num2str(batdate),'_',num2str(Bat),'.mat'));

% Load in spike data
load(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/B_ephys_data_aligned.mat'));

% Load in the audio event data (should include the wingbeats)
load(strcat(exp_data_path,'/audio/sorted_event_timestamps.mat'));
load(strcat(exp_data_path,'/audio/ttlConCat_segment_chunk_1_break_0.mat'));

min_pk_dist = 500000;
min_pk_ht = 0.35;

[pk_vals,pk_locs] = findpeaks(ttlConCat,'MinPeakDistance',min_pk_dist,'MinPeakHeight',min_pk_ht);
%pk_locs contains the indexes of the TTL pulses

% This is the value you will subtract from all vectors in audio space to
% zero it to the first TTL.
first_ttl_sample = pk_locs(1);
zerod_sorted_events = sorted_events-first_ttl_sample;

% First ephys TTL time
usec_of_first_ephys_TTL = (B_ephys_data.TT_unit(1).Timestamps(1)-B_ephys_data.TT_unit(1).AlignedTimestamps(1));
% ^^ This is the timestamp usec you want to subtract from all values in the CSC file to get the same alignment to the ephys spike file

% Load in the CSC.mat data
csc_files = dir(fullfile(strcat(exp_data_path,'ephys/logger',num2str(logger),'/extracted_data/*CSC*')));
load(strcat(csc_files(1).folder,'/',csc_files(1).name));
Fs = Estimated_channelFS_Transceiver(1);
raw_V = double(AD_count_int16*AD_count_to_uV_factor);
tstamps = Timestamps_of_first_samples_usec(1)+[0:length(raw_V)-1]/Fs*1e6;
tstamps_shift = tstamps-usec_of_first_ephys_TTL;   
t_stamps_g0 = find(tstamps_shift > 0); t_g0 = t_stamps_g0(1); t_stamps_g0=[];
tstamps = []; 
raw_V_shift = raw_V(t_g0:end); raw_V = [];

% % Plot the ciholas acceleration value for seg_amt samples ciholas seconds is /120
% % ephys seconds is /1e6
% seg_amt_sec = 780; % Take first 780 seconds of the session
% ciholas_sample_amt = seg_amt_sec*120;
% ephys_sample_amt = seg_amt_sec*1e6;

% Downsample the ephys data
t_sample_length = t(end)/length(t);
extra_time = tstamps_shift(end)/1e6-t(end);
extra_samples = round(extra_time/t_sample_length);
desired_ephys_trail_sample_length = length(t)+extra_samples;

oL = length(raw_V_shift);
y = interp1(1:oL,raw_V_shift,linspace(1,oL,desired_ephys_trail_sample_length)); y = y/3000;

% Plot ciholas data, wingbeats, and ephys!
Fs = 120;
figure(); hold on;
[up,lo] = envelope(a_flt+normrnd(0,1e-3,length(a_flt),1),Fs/10,'peak');     % Envelope of the acceleration signal (noise addedd to avoid problems with splining)
env = normalize(up - lo,'range');                                           % Amplitude of the envelope
env_th = otsuthresh(histcounts(env));                                       % Threshold (based on Otsu method). Can be set at 0.35
wBeats = movsum(env>env_th,1.67*Fs)>Fs/6; 
%figure(); hold on; %Heuristic criterion for flight detection
area(t,wBeats*3,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(t,normalize(v_abs,'range'));
plot(t,r(:,1),t,r(:,2));
plot(t,normalize(a_flt,'range',[-1 1]));
plot(t,normalize(movsum(env>env_th,1.67*Fs),'range'));
plot(t,up);
plot(t,lo);
plot(t,y(1:length(t)));
[y_up,y_lo] = envelope(y(1:length(t)),Fs/10,'peak');     % Envelope of the acceleration signal (noise addedd to avoid problems with splining)
y_env = normalize(y_up - y_lo,'range'); 
plot(t,y_up);
plot(t,y_lo);
hold off;

prompt = "Does this plot look aligned? (Y) or (N)";
g = input(prompt,'s');
if g=="Y"
    wingbeats_synced=1;
else
    wingbeats_synced=0;
end

% Calculate autocorreclaion between acceleration and ephys envelopes
start_xcorr_rng = [find(t==1159):find(t==1414)]; end_xcorr_rng = [find(t==5744):find(t==5825)];
[flight_ephys_xcorr,lags] = xcorr(y_up(start_xcorr_rng(1):start_xcorr_rng(end)),up(start_xcorr_rng(1):start_xcorr_rng(end)),10,'biased');
[flight_ephys_xcorr,lags] = xcorr(y(1:length(t)),a_flt);
figure(); stem(lags,flight_ephys_xcorr);

end