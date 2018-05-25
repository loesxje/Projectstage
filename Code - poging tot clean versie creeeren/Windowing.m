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
