clear, close, clc

reFs = 48000; %resampling rate in Hz
duration_window = 200; %in milliseconds

tic 

% Let op: De baseline wordt nu opgesteld aan de hand van de afwijkingsdata doordat dit voor nu even sneller werkt.
path = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\wavFiles\Dataset 1\afwijking_boor\'; 
filename_part_1 =  'mic_44100_s16le_channel_';
filename_part_2 = '_WAV.wav';
[baseline, feat_vec_total] = BaselineLib(reFs, duration_window, path, filename_part_1, filename_part_2);

toc

%%
clc

    subset = 5; % size of the subset
    part_mu_sigma = {};
    subset_mu_sigma = {};

for mic = 1:1%size(feat_vec_total,2) % over alle microfoons loopen
    feat_vec = feat_vec_total{mic};
    n_feat = size(feat_vec,1);
    for feat = 1:1%:n_feat                        % over alle features loopen
        one_feat = feat_vec(feat,:);
        n_windows = size(one_feat,2);
        for i = 1:subset:(n_windows - subset) % over alle featurevalues loopen
            subset_feat_vec = zeros(1, subset);
            for subset_i = 1:subset             % over de index in de subset loopen
                subset_feat_vec(subset_i) = one_feat(i + subset_i - 1);
            end
            subset_mu_sigma = {};
            mu = mean(subset_feat_vec);
            sigma = std(subset_feat_vec);
            subset_mu_sigma{feat, n_windows} = [mu; sigma]
        end
    end
end