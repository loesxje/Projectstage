clear, close, clc

reFs = 48000; %resampling rate in Hz
duration_window = 200; %in milliseconds
alpha = 1.5;

tic 

% Let op: De baseline wordt nu opgesteld aan de hand van de afwijkingsdata doordat dit voor nu even sneller werkt.
path = 'C:\Users\Loes\Documents\GitHub\Projectstage\wavFiles\Dataset 1\Normal\'; 
filename_part_1 =  'mic_44100_s16le_channel_';
filename_part_2 = '_WAV.wav';
[baseline, feat_vec_total] = BaselineLib(reFs, duration_window, path, filename_part_1, filename_part_2, alpha);

toc

%%
clc

    subset = 5; % size of the subset
    part_mu_sigma = {};
    subset_mu_sigma = {};

for mic = 1:1%size(feat_vec_total,2) % over alle microfoons loopen
    feat_vec = feat_vec_total{mic};
    n_feat = size(feat_vec,1);
        % over alle features loopen
        %one_feat = feat_vec(feat,:);
    n_windows = size(feat_vec,2);
    subset_feat_vec = {};%zeros(1, subset);
    n = 1;
    per_subset_feat = {};
    for i = 1:subset:(n_windows - subset) % neemt steeds subset van 5
        subset_feat_vec{n} = feat_vec(:,i:i+subset-1); % je hebt n sets, waarin rijen features zijn en kolommen bevatten steeds 5 windows
        n = n + 1; % uiteindelijk heb je n sets 
    end
    for j = 1:size(subset_feat_vec,2)
        per_subset = subset_feat_vec{j};
        per_subset_per_feat_mu_sigma = {};
        for k = 1:n_feat
            per_subset_per_feat = per_subset(k,:);
            mu = mean(per_subset_per_feat);
            sigma = std(per_subset_per_feat);
            per_subset_per_feat_mu_sigma{k} = [mu sigma];
        end
        per_subset_feat{j} = per_subset_per_feat_mu_sigma;
    end
end
