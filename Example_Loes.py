# -*- coding: utf-8 -*-
from Experiment_preprocessing import read_audio_data
from Experiment_preprocessing import stereo_to_mono
from Experiment_preprocessing import plot_raw_audio_data
from Experiment_preprocessing import split_vector

samplerate, data, times = read_audio_data('C:/Users/Loes/Documents/GitHub/Projectstage/audiovoorbeelden/heartbeat-noise-machine.wav')
mono = stereo_to_mono(data)
intervals = split_vector(times, mono)

# Fourier transform

from scipy.fftpack import dct

