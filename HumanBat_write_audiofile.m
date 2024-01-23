%% Writing an audiofile

filename_audio = 'feeder_click_ex2_bp.wav';
audioData = filteredAudio_2;
Fs = 192000;
audiowrite(filename_audio,audioData,Fs);

