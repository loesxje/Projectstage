%% Evaluate model using performance metrics
clc

% Run Anomaly_score for variable anom_score
% Run Anomaly_detection for variable detect_anom

%% Plot histogram of anomaly scores against their occurences per label
usual_scores = [];
unusual_scores = [];
for i = 1:numel(detect_anom) % Divide usual from unusual scores based on their labels
    if detect_anom(i) == 0
        usual_scores = [usual_scores, anom_score(i)];
    else
        unusual_scores = [unusual_scores, anom_score(i)];
    end
end

% plot histograms
histogram(usual_scores, 'NumBins', 100, 'FaceColor', 'b')
hold on
histogram(unusual_scores, 'NumBins', 100, 'FaceColor', 'r')

axis([0 100 0 100])
legend('usual', 'unusual')
xlabel('Anomaly score')
ylabel('Occurences')
%title('Performance Metrics')

%% Calculate RP@p
dist_prc = zeros(100,1);
for p = 1:100
    prc_usual = prctile(usual_scores, p);
    prc_unusual = prctile(unusual_scores,100-p);
    dist_prc(p) = vpa(prc_unusual - prc_usual);
end

%% Plot RP-curve

figure()
x = 1:1:100;
plot(x, dist_prc, 'Color', 'b', 'LineWidth', 1.25)
hold on
plot(x, zeros(1,100), 'Color', 'r')

xlabel('percentile')
ylabel('distance')
title('RP-Curve')
legend('Distance in percentiles')

%% Calculate RP-AUC 

a = 1;
b = 100;
m = ceil((b-a)/2);

Simpsons = ((b-a)/6) * (dist_prc(a) + 4*dist_prc(m) + dist_prc(b)) % Simpson's rule

txt = strcat('RP-AUC = ', int2str(Simpsons));
text(40, 75, txt) % Add RP-AUC value in plot



