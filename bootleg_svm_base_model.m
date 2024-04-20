function [SVMModel] = bootleg_svm_base_model()

    % Train an SVM model on pretty minimal data 

    % Load Table of values
    data = readtable('biosound_values.csv');
    
    % Create lable vector. 1 is echolocation, 0 is other.
    for i=1:length(data.calltype)
        if strcmp(data.calltype(i),'echo')
            y(i) = 1;
        elseif strcmp(data.calltype(i),'feeder')
            y(i) = 0;
        else
            y(i) = 2;
        end
    end
    
    % Insert labels into table.
    D = data(:,3:end);
    D.labels = y';
    
    % Make SVM
    features = table2array(D(:,[3]));
    labels = table2array(D(:,end));
    SVMModel = fitcecoc(features, labels);
    
    savedconfmat = {};
    num_runs = 10;
    for i = 1:num_runs
    
        % Cross-fold validation (partition the data into 1/5ths
        cv = cvpartition(size(features,1),'HoldOut',0.2);
        idx = cv.test;
        
        % Separate to training and test sets
        featuresTrain = features(~idx,:);
        labelsTrain = labels(~idx,:);
        featuresTest = features(idx,:);
        labelsTest = labels(idx,:);
        
        % Train the model on the training set (80% of data)
        SVMModel = fitcecoc(featuresTrain, labelsTrain);
        
        % Predict using the test set (20% of data)
        predictions = predict(SVMModel, featuresTest);
        
        % Evaluate the classifier (see how predictions match actual labels)
        confMat = confusionmat(labelsTest, predictions);
        accuracy = sum(diag(confMat)) / sum(confMat(:));
        fprintf('Accuracy: %.2f%%\n', accuracy * 100);
    
        savedconfmat{end+1} = confMat;
    end
end