clear, close, clc
%% Read audio file (baseline)
filename = 'C:\Users\Loes\Documents\GitHub\Projectstage\audiovoorbeelden\heartbeat-sounds\set_a\artifact__201105040918.wav';
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
% %determine duration of audio file
% duration = length(y)/reFs *1000; %(ms)
% duration_window = 200; %(ms)
% N = length(y); %array has N samples
% numWindows = floor(duration/duration_window); %number of windows you get
% N_perWindow = floor(N/numWindows); %how many samples will each window contains, door naar beneden af te ronden heb je niet precies 200 ms per frame meer
% 
% %Store samples in the corresponding window
% W = zeros(N_perWindow, numWindows); 
% begin = 1;
% stop =  N_perWindow;
% for i=1:numWindows
%     W(:,i) = y(begin:stop);
%     begin = stop+1;
%     stop = stop + N_perWindow;
% end
% %So W is the result of segmenting

%% Frequency spectrum extraction for each window
% W_freq = zeros(N_perWindow, numWindows);
% for j=1:numWindows
%     W_freq(:,j) = fft(W(:,j));
% end
% W_freq = abs(W_freq(2:N_perWindow/2+1,:)); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<aanpassing>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% %So W_freq is the result of frequency spectrum extracted from each window

%% Frequency spectrum extraction for whole signal
Freq = fftshift(y);
%% Calculate frequency bins with FFT <<<<<<<<<<<<<<<<<<<<<<<<< aanpassing
%Note that the index for the raw FFT are integers from 1?N.
%We need to process it to convert these integers to frequencies. That is where the sampling frequency counts. 
% N2 = N_perWindow;
N2 = length(y);
df = reFs/N2;
sampleIndex = 0:N2/2-1; %raw index for FFT plot
f = sampleIndex*df; %x-axis index converted to frequencies
%Now we can plot the absolute value of the FFT against frequencies as
subplot(2,1,1); stem(abs(Freq)); %sample values on x-axis
title('X[k]'); xlabel('k'); ylabel('|X(k)|');
subplot(2,1,2); plot(abs(Freq)); %x-axis represent frequencies, y-axis the magnitude response from FFT output
title('X[k]'); xlabel('frequencies (f)'); ylabel('|X(k)|');

%% Plot in time domain
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
% subplot(3,1,3)
% plot(f(1:Nfft/2),G(2:Nfft/2+1))


%% Feature extraction
feature_vector = FeaturesLib(y, Freq, f);

%% Reference vector Rescaled, i.e. normalized
maxRef = max(feature_vector);
minRef = min(feature_vector);
    
for e = 1:length(feature_vector)
    feature_vector(e) = (feature_vector(e) - minRef) / (maxRef - minRef);
end
 