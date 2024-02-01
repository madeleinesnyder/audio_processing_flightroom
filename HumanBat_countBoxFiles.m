%% Script to check how many of each kind of audio event we have

% Count how many of each kind of clip we have in the folder. 

audio_base = '/home/madeleine/Downloads/';
audio_files = dir(strcat(audio_base,'/*mic*'));

file_count = zeros(1,5);
for i=1:length(audio_files)
    audio_file_parts = split(audio_files(i).name,'_');
    if strcmp(audio_file_parts{1},'echolocation') | strcmp(audio_file_parts{1},'echo')
        file_count(1) = file_count(1)+1;
    elseif strcmp(audio_file_parts{1},'humanclick')
        file_count(2) = file_count(2)+1;
    elseif strcmp(audio_file_parts{1},'kevinvoice') | strcmp(audio_file_parts{1},'kevin')
        file_count(3) = file_count(3)+1;
    elseif strcmp(audio_file_parts{1},'madeleinevoice') | strcmp(audio_file_parts{1},'madeleine')
        file_count(4) = file_count(4)+1;
    elseif strcmp(audio_file_parts{1},'feederclick') | strcmp(audio_file_parts{1},'feeder')
        file_count(5) = file_count(5)+1;  
    end
end

disp(strcat(num2str(file_count(1))," echolocation snippits."));
disp(strcat(num2str(file_count(2))," humanclick snippits."));
disp(strcat(num2str(file_count(5))," feeder snippits."));
disp(strcat(num2str(file_count(3))," kevin voice snippits."));
disp(strcat(num2str(file_count(4))," madeleine voice snippits."));

