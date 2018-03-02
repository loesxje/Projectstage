
# -*- coding: utf-8 -*-
from Experiment_preprocessing import read_audio_data
from Experiment_preprocessing import plot_raw_audio_data
from Experiment_preprocessing import split_vector
from Experiment_preprocessing import fft_data

samplerate, mono, times = read_audio_data('C:/Users/Loes/Documents/GitHub/Projectstage/audiovoorbeelden/heartbeat-noise-machine.wav')
intervals = split_vector(times, mono)
Y_fft = fft_data(intervals)
    
   
from scipy.signal import spectrogram
samplerate, samples = wavfile.read('C:/Users/Loes/Documents/GitHub/Projectstage/audiovoorbeelden/heartbeat-noise-machine.wav')
freq, times, spectrogram = spectrogram(samples, samplerate)

plt.pcolormesh(times, freq, spectrogram)
plt.imshow(spectrogram)
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.show()

