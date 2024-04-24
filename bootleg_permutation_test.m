function [h_choice,h_choice_T] = bootleg_permutation_test(M_echo_hz,K_echo_hz,M_echo_hz_T,K_echo_hz_T)

    alpha = 0.05;
    % Perform a permutation test on the echolocation rate 
    disp("Performing permutation test")
    empirical_hz_diff = mean(M_echo_hz) - mean(K_echo_hz);
    empirical_hz_diff_T = mean(M_echo_hz_T) - mean(K_echo_hz_T);

    labeled_hz_matrix = [M_echo_hz,K_echo_hz;repmat(0,1,length(M_echo_hz)),repmat(1,1,length(K_echo_hz))]';
    labeled_hz_matrix_T = [M_echo_hz_T,K_echo_hz_T;repmat(0,1,length(M_echo_hz_T)),repmat(1,1,length(K_echo_hz_T))]';
    
    num_perms = factorial(length(labeled_hz_matrix(:,2)));
    if 1/num_perms > 0.02
        return 
    end

    shuffle_difference_values = [];
    shuffle_difference_values_T = [];
    temp = labeled_hz_matrix(:,2);
    for i=1:1000
        permute_indexes = randperm(length(temp));
        permute_array = temp(permute_indexes);
        values = labeled_hz_matrix(:,1);
        values_T = labeled_hz_matrix_T(:,1);
        s1 = values(permute_array==0);
        s2 = values(permute_array==1);
        diff_ = mean(s1)-mean(s2);
        shuffle_difference_values = [shuffle_difference_values,diff_];
        s1_T = values_T(permute_array==0);
        s2_T = values_T(permute_array==1);
        diff_T = mean(s1_T)-mean(s2_T);
        shuffle_difference_values_T = [shuffle_difference_values_T,diff_T];
    end

    figure(); hold on; title("Land: Null and empirial (red) for echolocation rate")
    histogram(shuffle_difference_values); xline(empirical_hz_diff,'r');
    figure(); hold on; title("Takeoff: Null and empirial (red) for echolocation rate")
    histogram(shuffle_difference_values_T); xline(empirical_hz_diff_T,'r');

    pval = sum(abs(shuffle_difference_values) >= abs(empirical_hz_diff)) / 1000;
    pval_T = sum(abs(shuffle_difference_values_T) >= abs(empirical_hz_diff_T)) / 1000;
    if pval >= alpha
        h_choice =  0;
    elseif pval < alpha
        h_choice = 1;
    end
    if pval_T >= alpha
        h_choice_T =  0;
    elseif pval_T < alpha
        h_choice_T = 1;
    end

end