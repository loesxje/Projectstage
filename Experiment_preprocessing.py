#-*- coding: utf-8 -*-
from scipy.io import wavfile
from matplotlib import pyplot as plt
import numpy as np
import math


# Load the data and calculate the time of each sample
samplerate, data = wavfile.read('C:/Users/Loes/Documents/GitHub/Projectstage/audiovoorbeelden/heartbeat-noise-machine.wav')

# Length of sampling Parameter
p = len(data)

times = np.arange(len(data[:p]))/float(samplerate)




# VECTORIZE DATA
# Stereo to Mono
def stereo_to_mono(audio_data):
    mono = []
    for i in range(len(audio_data[:p])):
        mono.append((audio_data[i,0] + data[i,1])/2);
        
    return mono

# =============================================================================
# # Make the plot
# # You can tweak the figsize (width, height) in inches
# plt.figure(figsize=(30, 4))
# #plt.fill_between(times, data[:p,0], data[:p,1], color='k') 
# plt.plot(times, mono[:p])
# plt.xlim(times[0], times[-1])
# plt.xlabel('time (s)')
# plt.ylabel('amplitude')
# # You can set the format by changing the extension
# # like .pdf, .svg, .eps
# plt.show()
# =============================================================================


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



# CALCULATE ENERGY



# (CREATE SPECTROGRAM)