%% Read audio file (baseline)
filename = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\Librosa\Chirping-Birds.wav'
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
N_perWindow = floor(N/numWindows) %how many samples will each window contains, door naar beneden af te ronden heb je niet precies 200 ms per frame meer

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

%So W_freq is the result of frequency spectrum extracted from each window

%% Plot in time domain
starttime = 0;
endingtime = duration/1000;
%t = 0:dt:(length(y)*dt)-dt;
t = linspace(starttime,endingtime,length(y)); %time vector
plot(t,y); xlabel('Seconds'); ylabel('Amplitude'); %Plot laat linker en rechter channel zien (blauw en rood)
%
figure
plot(psd(spectrum.periodogram,y,'Fs',reFs,'NFFT',length(y)));


%% plot in frequency domain
Nfft = length(y);
f = linspace(0,reFs,Nfft);
G = abs(fft(y(:,1), Nfft)); %the fft of the samples y in Nfft points
plot(f(1:Nfft/2),G(1:Nfft/2))


%% Feature extraction: Root mean square amplitude (time domain) 
W_squared = W.^2;
E_rms = zeros(1,numWindows);
for k = 1:numWindows
    E_rms(k) = sqrt(sum(W_squared(:,k))/length(W_squared));
end

%% Feature extraction: Amplitude variance (time domain)
E_var = zeros(1,numWindows);
for k = 1:numWindows
    window = W(:,k);
    mu = mean(window); %mean value of the window samples
    window_var = zeros(N_perWindow,1);
    for q = 1:N_perWindow
        window_var(q) = (window(q) - mu)^2;
    end
    window_var = sum(window_var);
    E_var(k) = window_var;
end

%% Feature extraction: Short time energy (time domain)
E_se = zeros(1,numWindows);
for i = 1:numWindows
    E_se(i) = sum(W(:,i).^2)
end
% so E_se is result of short time energy extracted from each window

%% Feature extraction: Zero Crossing Rate (time domain)
E_zcr = zeros(1,numWindows);
for i = 1:numWindows
    count = 0;
    for j = 1:N_perWindow
        if j < N_perWindow
            if W(j,i)>= 0 & W(j+1,i)< 0 || W(j,i)< 0 & W(j+1,i)>= 0
            count  = count + 1;
            end
        end
    end
    E_zcr(i) = count;
end

% so E_zcr is result of zero crossing rate extracted from each window

%% Feature extraction: Spectral peak (frequency domain)
E_sp = zeros(1,numWindows);
for i= 1:numWindows
    E_sp(i) = find(max(W_freq(:,i)));
end
