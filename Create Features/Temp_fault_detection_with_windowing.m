clc
%% Read audio file (baseline)
filename = 'C:\Users\Loes\Documents\GitHub\Projectstage\audiovoorbeelden\mic_44100_s16le_channel_8_WAV.wav';
[x,Fs] = audioread(filename); %The input values from audioread() are dimensionless, scaled to -1<=x<1
%sound(x,Fs) %play audio file
numChan = size(x,2); %number of channels where data comes from
x = x(:,1); %just pick one of the channel, assuming there is not much of time difference in signal
%% Resampling
reFs = 48000; %resampling rate
[p,q] = rat(reFs/Fs,0.0001); %find integers p and q that yield the correct resampling factor, tolerance of 0.0001 for resampling the signal that is very close to 48kHz
check = p/q*Fs; %check if the desired sample rate is obtained
y = resample(x,p,q); %y is array of resampled signal

%% Noise reduction
%% Segmenting
%determine duration of audio file
duration = length(y)/reFs *1000; %(ms)
duration_window = 200; %(ms)
N = length(y); %array has N samples
numWindows = floor(duration/duration_window); %number of windows you get
N_perWindow = floor(N/numWindows); %how many samples will each window contains, door naar beneden af te ronden heb je niet precies 200 ms per frame meer

%Store samples in the corresponding window
W = zeros(N_perWindow, numWindows); 
begin = 1;
stop =  N_perWindow;
for i=1:numWindows
    W(:,i) = y(begin:stop);
    begin = stop+1;
    stop = stop + N_perWindow;
end
%So W is the result of segmenting

%% Frequency spectrum extraction for each window
W_freq = zeros(N_perWindow, numWindows);
for j=1:numWindows
    W_freq(:,j) = fft(W(:,j));
end
W_freq = abs(W_freq(2:N_perWindow/2+1,:)); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<aanpassing>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%So W_freq is the result of frequency spectrum extracted from each window

%% Calculate frequency bins with FFT <<<<<<<<<<<<<<<<<<<<<<<<< aanpassing
%Note that the index for the raw FFT are integers from 1?N.
%We need to process it to convert these integers to frequencies. That is where the sampling frequency counts. 
N2 = N_perWindow;
df = reFs/N2;
sampleIndex = 0:N2/2-1; %raw index for FFT plot
f = sampleIndex*df; %x-axis index converted to frequencies
%Now we can plot the absolute value of the FFT against frequencies as
subplot(3,1,1); stem(sampleIndex,abs(W_freq(:,1))); %sample values on x-axis
title('X[k]'); xlabel('k'); ylabel('|X(k)|');
subplot(3,1,2); plot(f,abs(W_freq(:,1))); %x-axis represent frequencies, y-axis the magnitude response from FFT output
title('X[k]'); xlabel('frequencies (f)'); ylabel('|X(k)|');

%% Plot in time domain
%subplot(3,1,3)
starttime = 0;
endingtime = duration/1000;
%t = 0:dt:(length(y)*dt)-dt;
t = linspace(starttime,endingtime,length(y)); %time vector
plot(t,y); xlabel('Seconds'); ylabel('Amplitude'); %Plot laat linker en rechter channel zien (blauw en rood)
%
%figure
%plot(psd(spectrum.periodogram,y,'Fs',reFs,'NFFT',length(y)));


%% plot in frequency domain
% Nfft = length(y)-1;
% f = linspace(0,reFs,Nfft);
% G = abs(fft(y, Nfft)); %the fft of the samples y in Nfft points
% plot(f(1:Nfft/2),G(2:Nfft/2+1))

%% ALL IN ONE
clear, close, clc

reFs = 48000;
duration_window = 200;

[x, Fs] = ReadSignal('C:\Users\Loes\Documents\GitHub\Projectstage\audiovoorbeelden\', 'mic_44100_s16le_channel_8_WAV.wav');
[W, W_freq, numWindows, N_perWindow, df, f] = PreProcessingLibWindow(x, Fs, reFs, duration_window);
feature_vector_window = FeaturesLibWindow(W, numWindows, N_perWindow, W_freq, df, f);
    
%% Reference vector Rescaled, i.e. normalized
R_spread = mean(feature_vector_window(8,:));
R_flatness = mean(feature_vector_window(7,:));
R_centroid = mean(feature_vector_window(5,:));
R_sp = mean(feature_vector_window(5,:));
R_zcr = mean(feature_vector_window(4,:));
R_se = mean(feature_vector_window(3,:));
R_var = mean(feature_vector_window(2,:));
R_rms = mean(feature_vector_window(1,:));
Ref_vector  = [R_spread; R_flatness; R_centroid; R_sp; R_zcr; R_se; R_var; R_rms];
maxRef = max(Ref_vector);
minRef = min(Ref_vector);
    
for e = 1:length(Ref_vector)
    Ref_vector(e) = (Ref_vector(e) - minRef) / (maxRef - minRef);
end


