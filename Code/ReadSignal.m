function[x, Fs] = ReadSignal(path, filename)
    pathname = strcat(path, filename);
    [x,Fs] = audioread(pathname); %The input values from audioread() are dimensionless, scaled to -1<=x<1
    %sound(x,Fs) %play audio file
    numChan = size(x,2); %number of channels where data comes from
    x = x(:,1); %just pick one of the channel, assuming there is not much of time difference in signal
end