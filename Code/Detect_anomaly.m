clearvars -except baseline
close, clc

%load('baseline')
reFs = 48000; %resampling rate in Hz
duration_window = 200; %in milliseconds

% Afwijking path
path = 'C:\Users\Loes\Documents\GitHub\Projectstage\wavFiles\Dataset 1\afwijking_praten\'; 
filename_part_1 =  'mic_44100_s16le_channel_';
filename_part_2 = '_WAV.wav';
file = strcat(filename_part_1, int2str(1), filename_part_2); % hardcoded for file 0

% Afwijkingsdata inlezen, preprocessen en feature extraction uitvoeren
[x, Fs] = ReadSignal(path,file);
[W, W_freq, numWindows, N_perWindow, df, f] = PreProcessingLibWindow(x, Fs, reFs, duration_window);
feature_vector_window = FeaturesLibWindow(W, numWindows, N_perWindow, W_freq, df, f);

% deze f_v_w opdelen in dezelfde subsets als de baseline
    % baseline is opgedeeld per microfoon. Dus index 1 is voor mic 0, etc.
    % dus opdelen f_v_w in subsets is voor uitproberen nog niet nodig.
    % Later wel.

%% Stationaire standaarddeviatie baseline    
    
% de lb en ub uit de baseline trekken bij de juiste subset
bounds = baseline{1}(1,:);  % hardcoded voor feature 1 ( {1,:} )
lb = bounds(1);
ub = bounds(2);


% de f_v_w waarden met de lb en ub waarden vergelijken bij de juiste subset
result = zeros(length(feature_vector_window), length(baseline));    % lege vector om resultaat in op te slaan
anomaly = zeros(length(feature_vector_window), 1);

for i = 1:length(baseline)
    % de lb en ub uit de baseline trekken bij de juiste subset
    bounds = baseline{1}(i,:);  % {1} mic 1 -- {i,:} feature i
    lb = bounds(1);
    ub = bounds(2);
    for j = 1:length(feature_vector_window)
        % de f_v_w waarden eruit trekken 
        val = feature_vector_window(i,j); % f_v_w(i,j) met feature i en waarde j
        % de f_v_w waarden met de lb en ub waarden vergelijken bij de juiste subset
        if val > ub || val < lb
            result(j,i) = 1;
        end
    end
    for j = 1:length(feature_vector_window)
        % Raise anomaly if more than .. features are anomalous
        if sum(result(j,:) == 1) >= 3
            anomaly(j) = 1;
        end
    end
end


an = sum(anomaly == 1)
an/ length(feature_vector_window)
% Anomalie score berekenen

