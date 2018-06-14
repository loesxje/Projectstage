%Run eerst script semi_stationair.m
normaal = multi_feature_vectors;
MU1 = zeros(size(1,numFeatures));
SIGMA1 = zeros(numFeatures);
for f = 1:numFeatures
    MU1(f) = mean(normaal(f,:));
end
SIGMA1 = cov(multi_feature_vectors'); %%Dit moet de nieuwe SIGMAA zijn!!!
%% Run daarna eerst anomaly detection
afwijking = multi_feature_vectors_afwijking;

%pas op dat de naam hier zelfde is als bij 'normaal', houd daar rekening mee bij het runnen dus
p = vpa(mvnpdf(afwijking', MU1, SIGMA1))
%an anomaly als p<epsilon (zelf kiezen)

1/sqrt((2*pi)^8*norm(SIGMA1))*exp(-0.5*(
%zie voor interpretatie deze bron: https://nl.mathworks.com/help/stats/multivariate-normal-distribution.html
%uitleg over multivariate gaussian is hier te vinden: https://www.coursera.org/learn/machine-learning/lecture/DnNr9/anomaly-detection-using-the-multivariate-gaussian-distribution