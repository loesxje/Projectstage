
# -*- coding: utf-8 -*-
from Experiment_preprocessing import read_audio_data
from Experiment_preprocessing import plot_raw_audio_data
from Experiment_preprocessing import split_vector
from Experiment_preprocessing import fft_data

samplerate, data, times = read_audio_data('C:/Users/Loes/Documents/GitHub/Projectstage/audiovoorbeelden/heartbeat-sounds/set_a')
plot_raw_audio_data(times, data)
# =============================================================================
# intervals = split_vector(times, mono)
# Y_fft = fft_data(intervals)
#     
#    
# from scipy.signal import spectrogram
# samplerate, samples = wavfile.read('C:/Users/Loes/Documents/GitHub/Projectstage/audiovoorbeelden/heartbeat-noise-machine.wav')
# freq, times, spectrogram = spectrogram(samples, samplerate)
# 
# plt.pcolormesh(times, freq, spectrogram)
# plt.imshow(spectrogram)
# plt.ylabel('Frequency [Hz]')
# plt.xlabel('Time [sec]')
# plt.show()
# 
# =============================================================================
# =============================================================================
# from scipy.io import wavfile
# from matplotlib import pyplot as plt
# import numpy as np
# import math
# from numpy.fft import fft
# import os
# 
# folder = 'C:/Users/Loes/Documents/GitHub/Projectstage/audiovoorbeelden/heartbeat-sounds/set_a'
# 
# 
# files = [filename for filename in files[:5] for _, _, files in os.walk(folder)]
# samplerate, data = wavfile.read(os.path.join(folder, files[-4]) ) 
# times = np.arange(len(data[:]))/float(samplerate)
# data = data.astype(float)
# =============================================================================

