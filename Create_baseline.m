clear, close, clc
addpath('C:\Users\Loes\Documents\GitHub\Projectstage\Create Features')

reFs = 48000;
duration_window = 200;

tic 

% Let op. De baseline wordt nu opgesteld aan de hand van de afwijkingsdata Doordat dit voor nu even sneller werkt.
path = 'C:\Users\Loes\Documents\GitHub\Projectstage\wavFiles\Dataset 2\Afwijking_boor\'; 
filename_part_1 =  'mic_44100_s16le_channel_';
filename_part_2 = '_WAV.wav';
[baseline_cellarray, feat_vec_total] = BaselineLib(reFs, duration_window, path, filename_part_1, filename_part_2);

toc

%%
clc

    subset = 5; % size of the subset
    part_mu_sigma = {};
    subset_mu_sigma = {};
    t1 = 0;
    t2 = 0;
    
    for i = 1:1%numel(feat_vec_total) % over alle microfoons loopen
        feat_vec = feat_vec_total{i}; 
        subset_feat_vec = zeros(numel(feat_vec_total), subset);
        for j = 1:subset%numel(feat_vec(1,:)) % over 1 subset loopen
            one_feat = feat_vec(1,j);
            subset_feat_vec(1,j) = one_feat;
            t1 = t1 + 1;
            if t1 == subset
                range = zeros(2,size(subset_feat_vec,1));
                for k = 1:1%numel(feat_vec_total{i}) % over alle subsets loopen
                    range(1,k) = mean(subset_feat_vec(k,:));
                    range(2,k) = std(subset_feat_vec(k,:));
                    part_mu_sigma{i} = range
                end
            end
            % sizes goed krijgen. 
            %Daarna mu en sigma berekenen voor één
            % subset. 
            %Daarna voor meerdere subsets.
            %Dan voor alle features.
            % En tot slot voor alle mics.
            
%             range = zeros(2,8);
%             for j = 1:8
%                 range(1,j) = mean(feat_vec(j,:)); % mu (mean)
%                 range(2,j) = std(feat_vec(j,:)); % sigma (standarddeviation)
%                 part_mu_sigma{i} = range;
%             end
            t2 = t2 + 1;
            if t2 == 5
                subset_mu_sigma{i} = part_mu_sigma;
                t2 = 0;
            end
            part_mu_sigma = {};
        end
    end
    
%%

for mic = 1:size(feat_vec_total,2)
    feat_vec = feat_vec_total{mic};
    n_feat = size(feat_vec,1)
    subset_feat_vec = zeros(n_feat, subset);
    for feat = 1:n_feat
        one_feat = feat_vec(feat,:)
        for i = 1:size(one_feat,2)
            subset_feat_vec(feat) = one_feat(i)
        end
    end
end