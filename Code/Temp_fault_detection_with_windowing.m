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
feature_vector_window = FeaturesLibWindow(W, numWindows, N_perWindow, W_freq, df, f)
    
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


