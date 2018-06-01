clear, close, clc

reFs = 48000; %resampling rate in Hz
duration_window = 200; %in milliseconds

tic 

% Let op: De baseline wordt nu opgesteld aan de hand van de afwijkingsdata doordat dit voor nu even sneller werkt.
path = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\wavFiles\Dataset 1\Normal\'; 
filename_part_1 =  'mic_44100_s16le_channel_';
filename_part_2 = '_WAV.wav';
[baseline, feat_vec_total] = BaselineLib(reFs, duration_window, path, filename_part_1, filename_part_2);

toc

%%
clc
    subset = 5; % size of the subset, dus ongeveer 10 sec.
    part_mu_sigma = {};
    subset_mu_sigma = {};

for mic = 1:1%size(feat_vec_total,2) % over alle microfoons loopen
    feat_vec = feat_vec_total{mic};
    n_feat = size(feat_vec,1);
        %over alle features loopen
        %one_feat = feat_vec(feat,:);
    n_windows = size(feat_vec,2);
    subset_feat_vec = {};%zeros(1, subset);
    n = 1;
    per_subset_feat =[];
    for i = 1:subset:(n_windows - subset) % neemt steeds subset van 5
        subset_feat_vec{n} = feat_vec(:,i:i+subset-1); % je hebt n sets, waarin rijen features zijn en kolommen bevatten steeds 5 windows
        n = n + 1; % uiteindelijk heb je n sets 
    end
    for j = 1:size(subset_feat_vec,2)
        per_subset = subset_feat_vec{j};
        per_subset_per_feat_mu = zeros(8,1);
        for k = 1:n_feat
            per_subset_per_feat = per_subset(k,:);
            mu = mean(per_subset_per_feat);
            %sigma = std(per_subset_per_feat);
            per_subset_per_feat_mu(k) = mu;
        end
        per_subset_feat = [per_subset_feat per_subset_per_feat_mu];
    end
end

%dus nu heb je in per_subset_feat in de rijen de features, de kolommen
%sets. De waarde geeft de gemiddelde per feature per set aan 


%% Baseline opstellen
feat_mu = zeros(n_feat,1);
feat_sigma = zeros(n_feat,1);
for i = 1:n_feat
    feat_mu(i) = mean(per_subset_feat(i,:));
    feat_sigma(i) = std(per_subset_feat(i,:));
end
baseline_feat_mu_sigma = [feat_mu feat_sigma];
baseline_lower_upper = [baseline_feat_mu_sigma(:,1)-baseline_feat_mu_sigma(:,2) baseline_feat_mu_sigma(:,1)+baseline_feat_mu_sigma(:,2)]
%% afwijking inlezen

path = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\wavFiles\Dataset 1\afwijking_praten\'; 
filename_part_1 =  'mic_44100_s16le_channel_';
filename_part_2 = '_WAV.wav';
[baseline, feat_vec_total] = BaselineLib(reFs, duration_window, path, filename_part_1, filename_part_2);

% Dezelfde stap uitvoeren als hierboven: dus testset splitsen naar subsets
% van ... secondes
    subset = 5; % size of the subset, dus ongeveer 10 sec.
    part_mu_sigma = {};
    subset_mu_sigma = {};

for mic = 1:1%size(feat_vec_total,2) % over alle microfoons loopen
    feat_vec = feat_vec_total{mic};
    n_feat = size(feat_vec,1);
        %over alle features loopen
        %one_feat = feat_vec(feat,:);
    n_windows = size(feat_vec,2);
    subset_feat_vec = {};%zeros(1, subset);
    n = 1;
    per_subset_feat =[];
    for i = 1:subset:(n_windows - subset) % neemt steeds subset van 5
        subset_feat_vec{n} = feat_vec(:,i:i+subset-1); % je hebt n sets, waarin rijen features zijn en kolommen bevatten steeds 5 windows
        n = n + 1; % uiteindelijk heb je n sets 
    end
    for j = 1:size(subset_feat_vec,2)
        per_subset = subset_feat_vec{j};
        per_subset_per_feat_mu = zeros(8,1);
        for k = 1:n_feat
            per_subset_per_feat = per_subset(k,:);
            mu = mean(per_subset_per_feat);
            %sigma = std(per_subset_per_feat);
            per_subset_per_feat_mu(k) = mu;
        end
        per_subset_feat = [per_subset_feat per_subset_per_feat_mu];
    end
end

%%
afwijking_feat_mu = zeros(n_feat,1);
for i = 1:n_feat
    afwijking_feat_mu(i) = mean(per_subset_feat(i,:));
end
afwijking_feat_mu<baseline_lower_upper(1)  | afwijking_feat_mu>baseline_lower_upper(2)
    