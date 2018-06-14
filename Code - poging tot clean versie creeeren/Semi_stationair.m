%% Read audio data file
clear, clc
path = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\wavFiles\Dataset 1\normaal\'; 
filename = 'mic_44100_s16le_channel_0_TRAIN.wav';
[x, Fs] = ReadSignal(path, filename);

%% Resampling
reFs = 48000;
y = Resampling(reFs,Fs, x);

%% Windowing function helps split audio file in multiple input signals of 10 seconds each
duration_window = 60000; %ms
[S, numSamples, N_perSample] = Windowing(y, reFs, duration_window);

%% Then do the actual windowing for each input signal and extract features
duration_window = 200; %ms
multi_feature_vectors_normaal = [];
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
    multi_feature_vectors_normaal = [multi_feature_vectors_normaal feature_vector_signal];
end
 %% Normalize feature values for each featurevector
% for vect=1:numSamples
%     vector = multi_feature_vectors(:,vect);
%     maxRef = max(vector);
%     minRef = min(vector);
%     for elem = 1:numFeatures
%         vector(elem) = (vector(elem) - minRef) / (maxRef - minRef);
%     end
%     multi_feature_vectors(:,vect) = vector;
% end    

%% Normalize feature values for each featurevector according to x'i = (xi - mean)/std
% dit is niet meer goed, omdat het niet nodig is om de features te
% normaliseren
% for feat=1:numFeatures
%     meanFeat = mean(multi_feature_vectors(feat,:));
%     stdFeat = std(multi_feature_vectors(feat,:));
%     for instance=1:size(multi_feature_vectors,2);
%         multi_feature_vectors(feat,instance) = (multi_feature_vectors(feat,instance) - meanFeat)/stdFeat;
%     end
% end
%So multi_feature_vectors is now renewed with normalized featurevalues.
%% Calculate the Euclidic distance of each featurevector and store all
% multi_Eudist = zeros(1,numSamples);
% for vect = 1:numSamples
%    multi_Eudist(vect) = norm(multi_feature_vectors(:,vect));
% end
% 
% % Plot frequency histogram
% h = histogram(multi_Eudist);
% title('Frequency Histogram')
% ylabel('count')
% xlabel('Euclidic distance')
%h.NumBins = 10;

%% Create semi-stationair baseline based on Euclidic distance
% mean_Eudist = mean(multi_Eudist);
% std_Eudist = std(multi_Eudist);
% semi_stat_upper_lower = [mean_Eudist-std_Eudist mean_Eudist+std_Eudist];

%% Create semi-stationair baseline based on confidence interval
%  baseline_lower_upper = zeros(numFeatures,2);
%  for feat = 1:numFeatures
%      meanFeat = mean(multi_feature_vectors(feat,:));
%      stdFeat = std(multi_feature_vectors(feat,:));
%      baseline_lower_upper(feat,:) = [meanFeat-stdFeat meanFeat+stdFeat];
%  end

% %% Create semi-stationair baseline based on confidence interval (goede methode voor CI)
% baseline = zeros(numFeatures,2);
% mu_sigma = zeros(numFeatures,2);
% for feat = 1:numFeatures
%     meanFeat = mean(multi_feature_vectors_normaal(feat,:)); %mean of the samples
%     stdFeat = std(multi_feature_vectors_normaal(feat,:)); %standard deviation
%     mu_sigma(feat,:) = [meanFeat, stdFeat]; % save mu and sigma for every feat
%     n = length(multi_feature_vectors_normaal(feat,:)); %aantal samples
%     SEM = stdFeat/sqrt(n); %standard error
%     ts = tinv([0.025 0.975], n-1); %t-score
%     baseline(feat,:) = meanFeat + ts*SEM;
%     
% end
    
%% Store baseline in a .txt file
% fid = fopen('Baseline.txt','wt');
% for ii = 1:size(baseline,1)
%     fprintf(fid,'%20.18f \t',baseline(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid);

%% Run eerst script semi_stationar.m
MU = zeros(size(1,numFeatures));
SIGMA = zeros(numFeatures);
for f = 1:numFeatures
    MU(f) = mean(multi_feature_vectors_normaal(f,:));
    SIGMA(f,f) = std(multi_feature_vectors_normaal(f,:));
end
SIGMANIEUWWW = cov(multi_feature_vectors_normaal'); %%Dit moet de nieuwe SIGMAA zijn!!!