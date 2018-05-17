function[baseline, feat_vec_total] = BaselineLib(reFs, duration_window, path, filename_part_1, filename_part_2, alpha)
    feat_vec_total = create_total_feat_vec(reFs, duration_window, path, filename_part_1, filename_part_2);
    mu_sigma = calc_mu_sigma(feat_vec_total);
    baseline = calc_baseline(mu_sigma, feat_vec_total, alpha);
%   subset_mu_sigma = calc_subset_mu_sigma(feat_vec_total);
end

function[feat_vec_total] = create_total_feat_vec(reFs, duration_window, path, filename_part_1, filename_part_2)
   
    feat_vec_total = {};
    for i = 0:7
        nr = int2str(i);
        filename = strcat(filename_part_1, nr, filename_part_2);
        [x, Fs] = ReadSignal(path, filename);
        [W, W_freq, numWindows, N_perWindow, df, f] = PreProcessingLibWindow(x, Fs, reFs, duration_window);
        feature_vector_window = FeaturesLibWindow(W, numWindows, N_perWindow, W_freq, df, f);

        ind = i+1;
        feat_vec_total{ind} = feature_vector_window; % feat_vec_total{nr.feat, nr .wav}
    end
end

function[mu_sigma] = calc_mu_sigma(feat_vec_total)
    % Calculating the mu and sigma over the whole signal
    mu_sigma = {};

    for i = 1:numel(feat_vec_total)
        feat_vec = feat_vec_total{i};
        range = zeros(2,8);
        for j = 1:8
            range(1,j) = mean(feat_vec(j,:)); % mu (mean)
            range(2,j) = std(feat_vec(j,:)); % sigma (standarddeviation)    
        end
        mu_sigma{i} = range;
    end
end

function[baseline] = calc_baseline(mu_sigma, feat_vec_total, alpha)
    for i = 1:numel(feat_vec_total)
        mu_sigma{i}; % is van microfoon i de baseline voor alle features
        lower_upper_feat = zeros(8,2);
        for j = 1:numel(mu_sigma)
            mu = mu_sigma{i}(1,j); % is mu van microfoon i van j-de feature
            sigma = mu_sigma{i}(2,j); % is sigma van microfoon i van j-de feature
            lower_bound = mu - alpha * sigma; % is lower bound of the anomaly range
            upper_bound = mu + alpha * sigma; % is upper bound of the anomaly range
            lower_upper_feat(j,:) =  [lower_bound; upper_bound] %lower and upper bound per feature
        end
        baseline{i} = lower_upper_feat; % bounds with the different microfones as columns and the different features as rows
    end
end