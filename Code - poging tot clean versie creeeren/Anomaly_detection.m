%% Anomaly detection
% Read audio data file (same as baseline script but other path and file)
clc
path = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\wavFiles\Dataset 1\normaal\'; 
filename = 'mic_44100_s16le_channel_0_TEST.wav';
[x, Fs] = ReadSignal(path, filename);
%% Resampling (exact the same as baseline script)
reFs = 48000;
y = Resampling(reFs,Fs, x);

%% Windowing function helps split audio file in multiple input signals of 10 seconds each (exact the same as baseline script)
duration_window = 10000; %ms
[S, numSamples, N_perSample] = Windowing(y, reFs, duration_window);

%% Then do the actual windowing for each input signal and extract features (exact the same as baseline script)
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

%% Compare feature vector with the baseline (vervalt)
format long
comparison = zeros(size(multi_feature_vectors));

for s = 1:numSamples
    sample = multi_feature_vectors(:,s);
    for f = 1:numFeatures
        comparison(f,s) = sample(f)<baseline(f,1) | sample(f)>baseline(f,2);
    end
end
comparison       
        
%% Compare THE MEAN of the feature vectors with the baseline 
comparison = zeros(size(numFeatures,1));
    for f = 1:numFeatures
        gem = mean(multi_feature_vectors(f,:));
        comparison(f) = gem<baseline(f,1) | gem>baseline(f,2);
    end
comparison