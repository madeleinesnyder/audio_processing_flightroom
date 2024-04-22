function [h_choice] = bootleg_permutation_test(M_echo_hz,K_echo_hz)

    alpha = 0.05;
    % Perform a permutation test on the echolocation rate 
    disp("Performing permutation test")
    empirical_hz_diff = mean(M_echo_hz) - mean(K_echo_hz);

    labeled_hz_matrix = [M_echo_hz,K_echo_hz;repmat(0,1,length(M_echo_hz)),repmat(1,1,length(K_echo_hz))]';
    
    num_perms = factorial(length(labeled_hz_matrix(:,2)));
    if 1/num_perms > 0.02
        return 
    end

    shuffle_difference_values = [];
    temp = labeled_hz_matrix(:,2);
    for i=1:1000
        permute_indexes = randperm(length(temp));
        permute_array = temp(permute_indexes);
        values = labeled_hz_matrix(:,1);
        s1 = values(permute_array==0);
        s2 = values(permute_array==1);
        diff_ = mean(s1)-mean(s2);
        shuffle_difference_values = [shuffle_difference_values,diff_];
    end

    figure(); hold on;
    histogram(shuffle_difference_values); xline(empirical_hz_diff,'r');

    pval = sum(abs(shuffle_difference_values) >= abs(empirical_hz_diff)) / 1000;
    if pval >= alpha
        h_choice =  0;
    elseif pval < alpha
        h_choice = 1;
    end

end