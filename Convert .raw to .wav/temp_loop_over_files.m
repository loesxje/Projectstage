clear, close, clc
path = 'C:\Users\Loes\Documents\GitHub\Projectstage\Convert .raw to .wav\';
k = 7;
for i = 0:k
    variable = sprintf('mic_44100_s16le_channel_%d_WAV.wav',i)
    [x,Fs] = audioread(strcat(path,variable));
end