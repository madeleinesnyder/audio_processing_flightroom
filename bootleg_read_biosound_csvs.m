function [all_spectral_means,all_folder_timestamps] = bootleg_read_biosound_csvs(batdate)

% Go into box and make lists of the timestamps and the spectral mean of
% each observation

all_spectral_means = [];
all_folder_timestamps = [];
box_base = strcat('/home/',num2str(batdate),'/bootleg_snippits/');
box_base = strcat('/home/madeleine/Downloads/',num2str(batdate),'/bootleg_snippits/');
% For every folder in box base...
box_folders = dir(strcat(box_base,'/folder*'));
for i=1:length(box_folders)
    folder_path = strcat(box_base,box_folders(i).name,'/h5files/');
    wav_list = dir(strcat(folder_path,'/saved_event*'));
    biosound_csv_list = dir(strcat(folder_path,'/output*'));
    timestamps = [];
    for j=1:length(wav_list)
        temp_timestamp = split(wav_list(j).name,'_');
        timestamps = [timestamps;str2double(temp_timestamp{end}(1:end-3))];
    end
    try
        biosound_output = readtable(strcat(folder_path,biosound_csv_list(1).name));
    catch
        continue;
    end
    spectral_means = biosound_output.meanS;
    if size(biosound_output,1) ~= size(wav_list)
        disp("CRITICAL ERROR. Mismatch between samples in folder and samples in output file.");
        break;
    end

    all_folder_timestamps = [all_folder_timestamps;timestamps];
    all_spectral_means = [all_spectral_means;spectral_means];
end


