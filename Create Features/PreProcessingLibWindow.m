function[W_freq] = PreProcessingLibWindow(x, Fs, reFs, duration_window)
y = Resampling(reFs,Fs, x);    
[W, numWindows, N_perWindow] = Windowing(y, reFs, duration_window);
W_freq = FFT(W, numWindows, N_perWindow);
end


function[y] = Resampling(reFs, Fs, x)
    reFs = 48000; %resampling rate
    [p,q] = rat(reFs/Fs,0.0001); %find integers p and q that yield the correct resampling factor, tolerance of 0.0001 for resampling the signal that is very close to 48kHz
    check = p/q*Fs; %check if the desired sample rate is obtained
    y = resample(x,p,q); %y is array of resampled signal
end

function[W, numWindows, N_perWindow] = Windowing(y, reFs, duration_window)
    %determine duration of audio file
    duration = length(y)/reFs *1000; %(ms)
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
end

function[W_freq] = FFT(W, numWindows, N_perWindow)
    W_freq = zeros(N_perWindow, numWindows);
    for j=1:numWindows
        W_freq(:,j) = fft(W(:,j));
    end
    W_freq = abs(W_freq(2:N_perWindow/2+1,:)); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<aanpassing>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %So W_freq is the result of frequency spectrum extracted from each window
end

