clear, close, clc
path = 'C:\Users\Loes\Documents\GitHub\Projectstage\Convert .raw to .wav\';
k = 7;
i = 0;
variable = sprintf('mic_44100_s16le_channel_%d.wav',i)
[x,Fs] = audioread(strcat(path,variable))