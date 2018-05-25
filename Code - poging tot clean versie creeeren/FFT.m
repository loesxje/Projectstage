function[W_freq] = FFT(W, numWindows, N_perWindow)
    W_freq = zeros(N_perWindow, numWindows);
    for j=1:numWindows
        W_freq(:,j) = fft(W(:,j));
    end
    W_freq = abs(W_freq(2:N_perWindow/2+1,:));
    %So W_freq is the result of frequency spectrum extracted from each window
end