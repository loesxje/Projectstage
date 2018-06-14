%% Run daarna eerst anomaly detection
clc 

p = vpa(mvnpdf(multi_feature_vectors_afwijking', MU, SIGMA))
p = vpa(mvnpdf(multi_feature_vectors_normaal', MU, SIGMA))
%an anomaly als p<epsilon (zelf kiezen)




%zie voor interpretatie deze bron: https://nl.mathworks.com/help/stats/multivariate-normal-distribution.html
%uitleg over multivariate gaussian is hier te vinden: https://www.coursera.org/learn/machine-learning/lecture/DnNr9/anomaly-detection-using-the-multivariate-gaussian-distribution