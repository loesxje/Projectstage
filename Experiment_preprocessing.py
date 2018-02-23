# -*- coding: utf-8 -*-
import wave

# READ .WAV FILE
w = wave.open('C:/Users/Loes/Dropbox/Projectstage_Loes&Ngoc/audiovoorbeelden', 'r')
for i in range(w.getnframes()):
    frame = w.readframes(i)
    print frame


# VECTORIZE DATA



# SPLIT VECTOR



# FOURIER TRANSFORM DATA



# CALCULATE ENERGY



# (CREATE SPECTROGRAM)