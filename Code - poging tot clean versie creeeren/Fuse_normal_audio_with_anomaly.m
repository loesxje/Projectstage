%% Read audio: anomaly sound 1 minute
clear, clc
filename = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\wavFiles\Dataset 1\afwijking_boor\mic_44100_s16le_channel_0_WAV.wav';
[x,Fs] = audioread(filename); %The input values from audioread() are dimensionless, scaled to -1<=x<1
%sound(x,Fs) %play audio file
numChan = size(x,2); %number of channels where data comes from
x = x(:,1); %just pick one of the channel, assuming there is not much of time difference in signal
duration = length(x)/Fs *1000; %(ms)
%% Read audio: normal sound 10 minutes
filename = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\wavFiles\Dataset 1\normaal\mic_44100_s16le_channel_0_TRAIN.wav';
[y,Fs] = audioread(filename); %The input values from audioread() are dimensionless, scaled to -1<=x<1
%sound(x,Fs) %play audio file
numChan = size(x,2); %number of channels where data comes from
y = y(:,1); %just pick one of the channel, assuming there is not much of time difference in signal

%% Fuse audio
n = numel(x);
b = 44100*30; %30 sec
%m = 44100*5; %2 sec
normaal = y(n:n+n);
x(1:b) = normaal(1:b); %De eerste 30 seconden is omgevingsgeluid in normale staat en de overige 30 seconden is van afwijking
%x(b+m+1:end) = normaal(b+1:end);

%% Plot in time domain
starttime = 0;
endingtime = duration/1000;
%t = 0:dt:(length(y)*dt)-dt;
t = linspace(starttime,endingtime,length(x)); %time vector
plot(t,x); xlabel('Seconds'); ylabel('Amplitude'); %Plot laat linker en rechter channel zien (blauw en rood)
axis([-inf inf -20 20])

%% 
ax1 = subplot(211);
plot(x)
axis([-inf inf -10 10])
 
ax2 = subplot(212);
specgram(x)

linkaxes([ax1,ax2],'x')