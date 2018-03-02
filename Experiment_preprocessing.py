#-*- coding: utf-8 -*-
from scipy.io import wavfile
from matplotlib import pyplot as plt
import numpy as np
import math
from numpy.fft import fft
#from scipy.fftpack import dct


def read_audio_data(path):
    # Load the data and calculate the time of each sample
    samplerate, data = wavfile.read(path)    
    times = np.arange(len(data[:]))/float(samplerate)
    data = data.astype(float)
    # Turn audio_data from stereo to mono
    mono = data.sum(axis = 1) / 2
    
    return samplerate, mono, times


def plot_raw_audio_data(time, mono_data):
    # Make the plot
    # You can tweak the figsize (width, height) in inches
    plt.figure(figsize=(30, 4))
    #plt.fill_between(times, data[:p,0], data[:p,1], color='k') 
    plt.plot(time, mono_data[:])
    plt.xlim(time[0], time[-1])
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    # You can set the format by changing the extension
    # like .pdf, .svg, .eps
    plt.show()


# SPLIT VECTOR
def split_vector(time, mono_audio):
    interval_length = 0.02;
    n_splits = int(math.ceil(time[-1]/interval_length));
    interval_size = int(math.floor(len(mono_audio) / n_splits));
    
    intervals = []
    n = 0
    for i in range(n_splits+1):
        intervals.append(mono_audio[0+n:interval_size+n]);
        n += interval_size;
    
    return intervals


# FOURIER TRANSFORM DATA
    # bron: https://plot.ly/matplotlib/fft/
def fft_data(intervals):
    for i in range(len(intervals[0:10])):
        y = intervals[i+100];
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
      
# =============================================================================
#         fig, ax = plt.subplots(2,1)
#         ax[0].set_xlabel('Time')
#         ax[0].set_ylabel('Amplitude')
#         ax[0].plot(t,intervals[i+100])
#         ax[1].set_xlabel('Freq (Hz)')
#         ax[1].set_ylabel('|Y(freq)|')
#         ##### Not sure if E_F = F^2 or E_F = F
#         ax[1].plot(frq,abs(Y)**2, 'r')
# =============================================================================
    return Y

# CALCULATE ENERGY



# (CREATE SPECTROGRAM)
