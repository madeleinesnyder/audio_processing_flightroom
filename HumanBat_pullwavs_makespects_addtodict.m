%% Script to run through the Box folder and look at spectrograms. 
% Builds running dictionary of which filenames are good and uncorrupted to
% use

% Load in examples from Box
audio_base = '/home/madeleine/Downloads/';
audio_files = dir(strcat(audio_base,'/*mic*'));
for i=1:length(audio_files)
    audio_file_parts = split(audio_files(i).name,'_');
    if strcmp(audio_file_parts{1},'echolocation') | strcmp(audio_file_parts{1},'echo')
        sound_type = 'echo';
    elseif strcmp(audio_file_parts{1},'humanclick')
        sound_type = 'humanclick';
    elseif strcmp(audio_file_parts{1},'kevinvoice') | strcmp(audio_file_parts{1},'kevin')
        sound_type = 'kevin';
    elseif strcmp(audio_file_parts{1},'madeleinevoice') | strcmp(audio_file_parts{1},'madeleine')
        sound_type = 'madeleine';
    elseif strcmp(audio_file_parts{1},'feederclick') | strcmp(audio_file_parts{1},'feeder')
        sound_type = 'feeder';   
    end

    % Load in the file
    [data,fs] = audioread(strcat(audio_base,audio_files(i).name));
    HumanBat_basic_plots(data,fs,sound_type);

    user_qc_rating = [];
    userinput = input('Useable? "y" = yes, "n" = no :','s');
    if strcmp(userinput,'y')
        user_qc_rating = 1;
    elseif strcmp(userinput,'n')
        user_qc_rating = 0;
    else
        disp("Not a valid command. Excluding file");
        user_qc_rating = 0;
    end

    data_row = table({audio_files(i).name},{sound_type},[user_qc_rating],'VariableNames',headers);

    % If this looks good, add the filename ot a runing .mat directionary of
    % good sounds and save the dictionary back to box 
    if ~exist(strcat(audio_base,'audio_event_qc.csv'))
        disp("making audio_even_qc.csv file");
        %headers = {'filename','sound_type','qc'};
        %emptytable = array2table(cell(1,length(headers)));
        %emptytable.Properties.VariableNames = headers;
        csv_name = strcat(audio_base,'audio_event_qc.csv'); 
        %csv_data_updated = [data_row];
        writetable(data_row,csv_name);
    else
        csv_data = readtable(strcat(audio_base,'audio_event_qc.csv'),'Delimiter',',');
        csv_data_updated = [csv_data;data_row];
        newfilename = strcat(audio_base,'audio_event_qc.csv');
        writetable(csv_data_updated,newfilename);
    end
end
