#-*- coding: utf-8 -*-
import wave
import numpy as np

def _wav2array(nchannels, sampwidth, data):
    """data must be the string containing the bytes from the wav file."""
    num_samples, remainder = divmod(len(data), sampwidth * nchannels)
    if remainder > 0:
        raise ValueError('The length of data is not a multiple of '
                         'sampwidth * num_channels.')
    if sampwidth > 4:
        raise ValueError("sampwidth must not be greater than 4.")

    if sampwidth == 3:
        a = np.empty((num_samples, nchannels, 4), dtype=np.uint8)
        raw_bytes = np.fromstring(data, dtype=np.uint8)
        a[:, :, :sampwidth] = raw_bytes.reshape(-1, nchannels, sampwidth)
        a[:, :, sampwidth:] = (a[:, :, sampwidth - 1:sampwidth] >> 7) * 255
        result = a.view('<i4').reshape(a.shape[:-1])
    else:
        # 8 bit samples are stored as unsigned ints; others as signed ints.
        dt_char = 'u' if sampwidth == 1 else 'i'
        a = np.fromstring(data, dtype='<%s%d' % (dt_char, sampwidth))
        result = a.reshape(-1, nchannels)
    return result

 # READ .WAV FILE
w = wave.open('C:/Users/Loes/Documents/GitHub/Projectstage/audiovoorbeelden/074_heartbeat-noise-machine.wav', 'r')
rate = w.getframerate()
nchannels = w.getnchannels()
sampwidth = w.getsampwidth()
data = w.readframes(w.getnframes())
w.close()
array = _wav2array(nchannels, sampwidth, data)
print(array)
 
 




# VECTORIZE DATA



# SPLIT VECTOR



# FOURIER TRANSFORM DATA



# CALCULATE ENERGY



# (CREATE SPECTROGRAM)