%% Anomaly score berekenen
clc
    % Run Semi_stationair.m
    % Run Anomaly_detection.m
    % Run Multivariate_Gaussian.m


%% Calculate probability for anomaly score
probabilities_normaal = zeros(numFeatures, 1);
probabilities_afwijking = zeros(numFeatures, 1);

for feat = 1:numFeatures
    %pd = makedist('Normal', mu_sigma(feat,1), mu_sigma(feat,2));
    val = normaal(feat);
    probabilities_normaal(feat) = mvnpdf(val, mu_sigma(feat,1), mu_sigma(feat,2)); % calculate probability density function of normal dist %icdf(pd, val);
    val = afwijking(feat);
     probabilities_afwijking(feat) = mvnpdf(val, mu_sigma(feat,1), mu_sigma(feat,2)); % calculate probability density function of normal dist %icdf(pd, val);
end
probabilities_normaal;
probabilities_afwijking

%% Normalize probabilities
minimum = min(probabilities_normaal(:)); % min of P of the feat of the anomaly data
maximum = max(probabilities_normaal(:)); % max of P of the feat of the anomaly data
normalized = zeros(numFeatures,1);
for feat = 1:numFeatures
    normalized(feat) = 1 - (probabilities_afwijking(feat) - minimum) / (maximum - minimum) ; % min-max normalization
end
anom_score = normalized * 100

% De min/max is alleen genomen over de anomaly data en niet ook over de normale data
% Verder vraag ik mij nog af waarom sommige waarden bij probabilities
% groter zijn dan 1. Ik dacht dat p-values niet groter konden zijn dan 1.
% Zoek dit nog even uit.

%   MIN EN MAX VALUES ZIJN ALLEEN GEBASEERD OP DE ANOMALY DATA. BASEER ZE
%   OOK OP DE NORMALE DATA.

%% Probeerscript
% x = [-2,-1,0,1,2,3,4,5,6];
% mu = 2;
% sigma = 1;
% y = normpdf(x,mu,sigma)'
% 
% minimum = min(y);
% maximum = max(y);
% 
% normalized = zeros(length(y),1);
% for feat = 1:length(y)
%     normalized(feat) = 1 - (y(feat) - minimum) / (maximum - minimum); % min-max normalization
% end
% normalized