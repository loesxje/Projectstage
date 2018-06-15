%% Anomaly detection
% Read audio data file (same as baseline script but other path and file)
            %clc
path = 'C:\Users\Loes\Documents\GitHub\Projectstage\wavFiles\Dataset 1\normaal\'; 
filename = 'mic_44100_s16le_channel_0_TEST.wav';
[x, Fs] = ReadSignal(path, filename);
%% Resampling (exact the same as baseline script)
reFs = 48000;
y = Resampling(reFs,Fs, x);

%% Windowing function helps split audio file in multiple input signals of 10 seconds each (exact the same as baseline script)
duration_window = 5000; %ms
[S, numSamples, N_perSample] = Windowing(y, reFs, duration_window);

%% Then do the actual windowing for each input signal and extract features (exact the same as baseline script)
duration_window = 200; %ms
multi_feature_vectors_afwijking = [];
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
    multi_feature_vectors_afwijking = [multi_feature_vectors_afwijking feature_vector_signal];
end


%% Compare data to baseline (detect anomaly)
probabilities_afwijking = vpa(mvnpdf(multi_feature_vectors_afwijking', MU, VAR)); % pdf value

alpha = 3; % How many STD's from the mean
p_bound = vpa(mvnpdf((MU' - (alpha * diag(STD)))', MU, VAR)); % p value for the lower and upper bounds
detect_anom = [];

for p = 1:numel(probabilities_afwijking)
    if (probabilities_afwijking(p) < double(p_bound))
        detect_anom = [detect_anom 1];
    else
        detect_anom = [detect_anom 0];
    end
end
detect_anom

