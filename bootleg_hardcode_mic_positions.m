function [mic_positions] = bootleg_hardcode_mic_positions()

    batdate = 221125;

    potential_mic_coords = [-2800 0 1500; 0 2600 1500; 2800 0 1500; 0 -2600 1500];

    exp_data_path = strcat('/home/madeleine/mnt/server2/users/KQMS/HumanBat/1464314684/processed/',num2str(batdate),'/');

    % load in human data 
    load(strcat(exp_data_path,'ciholas/aligned_human_position_data','.mat'));
    madeleine = 2; kevin = 4;

    % at a certain time plot where humans are (1hr and 2minutes = )
    timestamp_in_sec = 62*60;
    timestamp_in_ciholas = timestamp_in_sec*120;

    % Plot positions 
    figure(); hold on; xlim([-2900 2900]); ylim([-2600 2600]); zlim([0 2300]);
    plot3(human_r(timestamp_in_ciholas-120:timestamp_in_ciholas+120*10,1,madeleine),human_r(timestamp_in_ciholas-120:timestamp_in_ciholas+120*10,2,madeleine),human_r(timestamp_in_ciholas-120:timestamp_in_ciholas+120*10,3,madeleine),'g');
    plot3(human_r(timestamp_in_ciholas-120:timestamp_in_ciholas+120*10,1,kevin),human_r(timestamp_in_ciholas-120:timestamp_in_ciholas+120*10,2,kevin),human_r(timestamp_in_ciholas-120:timestamp_in_ciholas+120*10,3,kevin),'r');
    text(potential_mic_coords(1,1),potential_mic_coords(1,2), ['mic2'],'FontWeight', 'bold');
    text(potential_mic_coords(2,1),potential_mic_coords(2,2), ['mic3'],'FontWeight', 'bold');
    text(potential_mic_coords(3,1),potential_mic_coords(3,2), ['mic4'],'FontWeight', 'bold');
    text(potential_mic_coords(4,1),potential_mic_coords(4,2), ['mic1'],'FontWeight', 'bold');

    mic_positions = [potential_mic_coords(4,:); potential_mic_coords(1,:); potential_mic_coords(2,:); potential_mic_coords(3,:)];
end