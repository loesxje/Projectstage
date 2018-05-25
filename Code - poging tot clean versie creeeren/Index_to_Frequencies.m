function[df, f] = Index_to_Frequencies(reFs, N_perWindow)
    df = reFs/N_perWindow;
    sampleIndex = 0:N_perWindow/2-1; %raw index for FFT plot
    f = sampleIndex*df; %x-axis index converted to frequencies
end