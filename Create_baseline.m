clear, close, clc
addpath('C:\Users\Loes\Documents\GitHub\Projectstage\Create Features')

reFs = 48000;
duration_window = 200;

feat_vec_total = {};
for i = 0:7
    nr = int2str(i);
    filename = strcat( 'mic_44100_s16le_channel_', nr, '_WAV.wav');
    [x, Fs] = ReadSignal('C:\Users\Loes\Documents\GitHub\Projectstage\wavFiles\Waterkoker\', filename);
    [W, W_freq, numWindows, N_perWindow, df, f] = PreProcessingLibWindow(x, Fs, reFs, duration_window);
    feature_vector_window = FeaturesLibWindow(W, numWindows, N_perWindow, W_freq, df, f);
            
    ind = i+1;
    feat_vec_total{ind} = feature_vector_window; % feat_vec_total{1, nr .wav}
end


%% make baseline
% bepaal mu en sigma per feature om zo een range (baseline) te
% determineren

mu_sigma = {};
baseline = {};

for i = 1:numel(feat_vec_total)
    feat_vec = feat_vec_total{i};
    range = zeros(2,8);
    for j = 1:8
        range(1,j) = mean(feat_vec(j,:)); % mu (mean)
        range(2,j) = std(feat_vec(j,:)); % sigma (standarddeviation)
        mu_sigma{i} = range;
    end
end

for i = 1:numel(feat_vec_total)
    mu_sigma{i}; % is van microfoon i de baseline voor alle features
    for j = 1:numel(mu_sigma)
        mu = mu_sigma{i}(1,j); % is mu van microfoon i van j-de feature
        sigma = mu_sigma{i}(2,j); % is sigma van microfoon i van j-de feature
        alpha = 1;
        lower_bound = mu - alpha * sigma; % is lower bound of the anomaly range
        upper_bound = mu + alpha * sigma; % is upper bound of the anomaly range
        baseline{j, i} = [lower_bound; upper_bound]; % bounds with the different microfones as colomns and the different features as rows.
    end
 
end
