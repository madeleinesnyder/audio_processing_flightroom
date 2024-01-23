function [y] = hbu_wienerFilter(x, w)
% Weiner filter where x is the signal and w is the filter window size
    N = length(x);
    y = zeros(size(x));
    for i = 1:N
        if i == round(length(x)/6)
            disp("sixth o the way there")
        end
        idx = max(1, i-w):min(N, i+w);
        localMean = mean(x(idx));
        localVariance = var(x(idx));
        noiseVariance = 0.01^2; % Change according to noise level, if known
        y(i) = localMean + (max(0, localVariance - noiseVariance) / localVariance) * (x(i) - localMean);
    end
end