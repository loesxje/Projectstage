#-*- coding: utf-8 -*-
import wave
import struct
# READ .WAV FILE
# =============================================================================
# w = wave.open('C:/Users/Loes/Documents/GitHub/Projectstage/audiovoorbeelden/074_heartbeat-noise-machine.wav', 'r')
# for i in range(w.getnframes()):
#     print i
#     frame = w.readframes(i)
# =============================================================================
aud = open("C:/Users/Loes/Documents/GitHub/Projectstage/audiovoorbeelden/074_heartbeat-noise-machine.wav", "r")

data = []

for i in range(len(aud.read())):
    SnString = aud.read(2)
    Sn = struct.unpack("h", SnString)
    data.append(Sn)
print data    


# VECTORIZE DATA



# SPLIT VECTOR



# FOURIER TRANSFORM DATA



# CALCULATE ENERGY



# (CREATE SPECTROGRAM)