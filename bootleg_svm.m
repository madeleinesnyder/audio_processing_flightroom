function [valid_echolocation_timestamps] = bootleg_svm(spectral_means,SVMModel,timestamps)

    % Returns whether each sample is an echolocation or not. 
    disp("Getting echo predictions")
    predictions = predict(SVMModel, spectral_means);
    valid_echolocation_timestamps = timestamps(predictions==0);

end
