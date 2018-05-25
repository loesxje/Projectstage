function[y] = Resampling(reFs, Fs, x)
    %reFs = 48000; %resampling rate ->>kan weg vgm
    [p,q] = rat(reFs/Fs,0.0001); %find integers p and q that yield the correct resampling factor, tolerance of 0.0001 for resampling the signal that is very close to 48kHz
    check = p/q*Fs; %check if the desired sample rate is obtained
    y = resample(x,p,q); %y is array of resampled signal
end
