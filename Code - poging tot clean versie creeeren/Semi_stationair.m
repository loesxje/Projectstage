%% Read audio data file
clear, clc
path = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\wavFiles\Dataset 1\Normal\'; 
filename = 'mic_44100_s16le_channel_0_TRAIN.wav';
[x, Fs] = ReadSignal(path, filename);

%% Resampling
reFs = 48000;
y = Resampling(reFs,Fs, x);

%% Windowing function helps split audio file in multiple input signals of 10 seconds each
duration_window = 10000; %s
[S, numSamples, N_perSample] = Windowing(y, reFs, duration_window);

%% Then do the actual windowing for each input signal and extract features
duration_window = 200; %ms
multi_feature_vectors = [];
for s=1:numSamples %for each input signal
    
    %Window the signal to 200 ms consecutive frames 
    y = S(:,s);
    [W, numWindows, N_perWindow] = Windowing(y, reFs, duration_window);
    
    %Transform each window to fft
    W_freq = FFT(W, numWindows, N_perWindow);
    
    % Index to frequency
    [df, f] = Index_to_Frequencies(reFs, N_perWindow);
    
    %Extract features from each window
    features_window = FeaturesLibWindow(W, numWindows, N_perWindow, W_freq, df, f);
    
    %Calculate mean feature values for each input signal: obtaining a featurevector of the input signal
    numFeatures = size(features_window,1);
    feature_vector_signal = zeros(numFeatures,1);
    for feat=1:numFeatures
        feature_vector_signal(feat) = mean(features_window(feat,:));
    end
    
    %Store all featurevectors of numSamples
    multi_feature_vectors = [multi_feature_vectors feature_vector_signal];
end
%% Normalize feature values for each featurevector
for vect=1:numSamples
    vector = multi_feature_vectors(:,vect);
    maxRef = max(vector);
    minRef = min(vector);
    for elem = 1:numFeatures
        vector(elem) = (vector(elem) - minRef) / (maxRef - minRef);
    end
    multi_feature_vectors(:,vect) = vector;
end    

%% Calculate the Euclidic distance of each featurevector and store all
multi_Eudist = zeros(1,numSamples);
for vect = 1:numSamples
   multi_Eudist(vect) = norm(multi_feature_vectors(:,vect));
end

% Plot frequency histogram
h = histogram(multi_Eudist);
title('Frequency Histogram')
ylabel('count')
xlabel('Euclidic distance')
%h.NumBins = 10;

%% Create semi-stationair baseline
mean_Eudist = mean(multi_Eudist);
std_Eudist = std(multi_Eudist);
semi_stat_upper_lower = [mean_Eudist-std_Eudist mean_Eudist+std_Eudist];

