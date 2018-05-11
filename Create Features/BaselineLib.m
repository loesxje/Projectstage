function[baseline, feat_vec_total] = BaselineLib(reFs, duration_window, path, filename_part_1, filename_part_2)
    feat_vec_total = create_total_feat_vec(reFs, duration_window, path, filename_part_1, filename_part_2);
    mu_sigma = calc_mu_sigma(feat_vec_total);
    baseline = calc_baseline(mu_sigma, feat_vec_total);
%     subset_mu_sigma = calc_subset_mu_sigma(feat_vec_total);
end

function[feat_vec_total] = create_total_feat_vec(reFs, duration_window, path, filename_part_1, filename_part_2)
    addpath('C:\Users\Loes\Documents\GitHub\Projectstage\Create Features')
    
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
            mu_sigma{i} = range;
        end
    end
end

function[baseline] = calc_baseline(mu_sigma, feat_vec_total)
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
end

% function[subset_mu_sigma] = calc_subset_mu_sigma(feat_vec_total)
%     part_mu_sigma = {};
%     subset_mu_sigma = {};
%     t = 0;
%     
%     for i = 1:numel(feat_vec_total)
%         feat_vec = feat_vec_total{i};
%         for k = 1:numel(feat_vec{i}(;,1))
%             
%             
%             range = zeros(2,8);
%             for j = 1:8
%                 range(1,j) = mean(feat_vec(j,:)); % mu (mean)
%                 range(2,j) = std(feat_vec(j,:)); % sigma (standarddeviation)
%                 part_mu_sigma{i} = range;
%             end
%             t = t + 1;
%             if t == 5
%                 subset_mu_sigma{i} = part_mu_sigma;
%                 t = 0;
%             end
%             part_mu_sigma = {};
%         end
%     end
% end