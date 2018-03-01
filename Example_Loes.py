
# -*- coding: utf-8 -*-
from Experiment_preprocessing import read_audio_data
from Experiment_preprocessing import stereo_to_mono
from Experiment_preprocessing import plot_raw_audio_data
from Experiment_preprocessing import split_vector
import numpy as np
import matplotlib.pyplot as plt

samplerate, data, times = read_audio_data('C:/Users/Loes/Documents/GitHub/Projectstage/audiovoorbeelden/heartbeat-noise-machine.wav')
mono = stereo_to_mono(data)
intervals = split_vector(times, mono)

# Fourier transform

from numpy.fft import fft

for i in range(len(intervals[0:5x])):
    y = intervals[i];
    Fs = len(y); # sampling rate
    Ts = 1.0/Fs; # sampling interval
    t = np.arange(0,1,Ts) # time vector
    ff = 5; # freq of the signal
    n = len(y);
    k = np.arange(n);
    T = n/Fs
    frq = k/T; # two sides freq range
    frq = frq[range(n/2)] # one side freq
    Y = fft(y); # fft computing and normalization
    Y = Y[range(n/2)];
    
    fig, ax = plt.subplots(2,1)
    ax[0].set_xlabel('Time')
    ax[0].set_ylabel('Amplitude')
    ax[0].plot(t,y)
    ax[1].set_xlabel('Freq (Hz)')
    ax[1].set_ylabel('|Y(freq)|')
    ax[1].plot(frq,abs(Y), 'r')
    
    


