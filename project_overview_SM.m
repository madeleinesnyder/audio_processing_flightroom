%% Project

%% 1. Gather .wav data for training data

% DONE (MAREESA)

%% 2. Get bioSound values for training data

% DONE (SANJANA)

%% 3. Train SVM on hand selected data

[SVMModel] = bootleg_svm_base_model();


%% 4. BATCH PROCESSING!! (B-steps):

    %% B1. Extract all "events" from audio signal in janky way
    
    bootleg_findpeaks.m 
    
    %% B2. Turn these timestamps into wav files (upload to box)
    
    bootleg_make_wavs(Bat,batdate);
    
    %% B3 Run these through biosound and get the features for each sample (Biosound)
    
    %% B4. Read out the biosound measures for each of these audio events 
    
    [spectral_means,timestamps] = bootleg_read_biosound_csvs(batdate);
    
    %% B5. Use SVM to classify this batch of samples as yes or no echolocation. Get valid timestamps. SAVE.
    
    predictions = predict(SVMModel, spectral_means);
    valid_echolocation_timestamps = timestamps(predictions==1);
    valid_timestamp_filename = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/audio/saved_echolocation_timestamps_biosound_verified.mat');
    save(valid_timestamp_filename,'valid_echolocation_timestamps');

    %% B6. Look at random samples to see if they make sense. 
    %bootleg_plot_sample(batdate,timestamps)



   
%% PLOTTING WITH FLIGHT DATA

%% 5. Plot the echolocations on the flight trajectories (start buffer is in samples)
Bat = 32626; batdate=221128; logger=15; start_buffer=0; end_buffer=0; cluster = 4; plot_units = 1; unit=1;
[M_flight_hz_samples_L,K_flight_hz_samples_L,M_flight_hz_samples_T,K_flight_hz_samples_T] = bootleg_plot_echos_on_flights(Bat,batdate,logger,unit,start_buffer,end_buffer,cluster,plot_units)

%% Perform permutation test to see if the echolocation rate at landing is significantly different between flights to M and flights to K
[sig_diff_L,sig_diff_T] = bootleg_permutation_test(M_flight_hz_samples_L,K_flight_hz_samples_L,M_flight_hz_samples_T,K_flight_hz_samples_T);

if sig_diff_L == 1
    disp("Significnat difference in mean number of echolocations upon landing to M v.s. to K!");
else
    disp("No sig diff in echolocation rate upon landing on M v.s. K");
end


if sig_diff_T == 1
    disp("Significnat difference in mean number of echolocations upon takeoff to M v.s. to K!");
else
    disp("No sig diff in echolocation rate upon takoff on M v.s. K");
end

