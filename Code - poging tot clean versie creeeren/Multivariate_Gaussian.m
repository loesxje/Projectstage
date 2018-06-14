%Run eerst script semi_stationar.m
normaal = multi_feature_vectors;
MU1 = zeros(size(1,numFeatures));
SIGMA1 = zeros(numFeatures);
for f = 1:numFeatures
    MU1(f) = mean(normaal(f,:));
    SIGMA1(f,f) = std(normaal(f,:));
end
%% Run daarna eerst anomaly detection
afwijking = multi_feature_vectors;

%pas op dat de naam hier zelfde is als bij 'normaal', houd daar rekening mee bij het runnen dus
p = vpa(mvnpdf(afwijking', MU1, SIGMA1))
%an anomaly als p<epsilon (zelf kiezen)


%zie voor interpretatie deze bron: https://nl.mathworks.com/help/stats/multivariate-normal-distribution.html
%uitleg over multivariate gaussian is hier te vinden: https://www.coursera.org/learn/machine-learning/lecture/DnNr9/anomaly-detection-using-the-multivariate-gaussian-distribution