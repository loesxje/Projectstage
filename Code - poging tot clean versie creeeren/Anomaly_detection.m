%% Anomaly detection
% Read audio data file (same as baseline script but other path and file)
            %clc
path = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\wavFiles\Dataset 1\afwijking_praten\'; 
filename = 'mic_44100_s16le_channel_0_WAV.wav';
[x, Fs] = ReadSignal(path, filename);
%% Resampling (exact the same as baseline script)
reFs = 48000;
y = Resampling(reFs,Fs, x);

%% Windowing function helps split audio file in multiple input signals of 10 seconds each (exact the same as baseline script)
duration_window = 10000; %ms
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
tic
probabilities_afwijking = mvnpdf(multi_feature_vectors_afwijking', MU, VAR); % pdf value

alpha = 3; % How many STD's from the mean
p_bound = mvnpdf((MU' - (alpha * diag(STD)))', MU, VAR); % p value for the lower and upper bounds
detect_anom = [];

for p = 1:numel(probabilities_afwijking)
    if probabilities_afwijking(p) < p_bound
        detect_anom = [detect_anom 1];
    else
        detect_anom = [detect_anom 0];
    end
end
detect_anom

toc

%% Original Gaussian model (IS NIET SNELLER) 
tic
px = zeros(size(multi_feature_vectors_afwijking, 2),1); %p(x) value for each samples to judge whether it is an anomaly or not
for k = 1:size(multi_feature_vectors_afwijking,2)
    px_feat = zeros(1,numFeatures);
    for f = 1:numFeatures
        px_feat(f) = normpdf(multi_feature_vectors_afwijking(f,k), MU(f), STD(f));
    end
    px(k) = prod(px_feat);
end        
detect_anomaly = [];
for p = 1:size(multi_feature_vectors_afwijking,2)
    if px(p) < p_bound
        detect_anomaly = [detect_anomaly 1];
    else
        detect_anomaly = [detect_anomaly 0];
    end
end
detect_anomaly 
toc