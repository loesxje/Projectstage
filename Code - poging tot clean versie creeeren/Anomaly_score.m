%% Anomaly score berekenen
clc
    % Run Semi_stationair.m
    % Run Anomaly_detection.m
    

%% Calculate probability for anomaly score wiht Multivariate Gaussian
% If the values are above the MU, then calculate the 1 - cdf(val, mu, sigma)
% If below the MU, then calculate cdf(val, mu, sigma)

% for i = multi_feature_vectors_afwijking
%     if mean(i) > mean(MU)
%         probabilities_afwijking = 1 - vpa(mvncdf(multi_feature_vectors_afwijking', MU, SIGMA));
%     else
%         probabilities_afwijking = vpa(mvncdf(multi_feature_vectors_afwijking', MU, SIGMA));
%     end
% end
% 
% for i = multi_feature_vectors_normaal
%     if mean(i) > mean(MU)
%         probabilities_normaal = 1 - vpa(mvncdf(multi_feature_vectors_normaal', MU, SIGMA));
%     else
%         probabilities_normaal = vpa(mvncdf(multi_feature_vectors_normaal', MU, SIGMA));
%     end
% end

probabilities_afwijking = vpa(mvncdf(multi_feature_vectors_afwijking', MU, VAR)); % cdf value
probabilities_normaal = vpa(mvncdf(multi_feature_vectors_normaal', MU, VAR)); % cdf value


%% Normalize probabilities
minimum = min([min(probabilities_afwijking), min(probabilities_normaal)]); % min of multivariate gaussian distribution of the normal data
maximum = max([max(probabilities_afwijking), max(probabilities_normaal)]); % max of multivariate gaussian distribution of the normal data
normalized = zeros(length(probabilities_afwijking),1);
for i = 1:length(probabilities_afwijking)
    normalized(i) = (probabilities_afwijking(i) - minimum) / (maximum - minimum) ; % min-max normalization
end

%% Calculate anomaly score
anom_score = normalized * 100

