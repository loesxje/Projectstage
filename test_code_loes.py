import lib_Loes as ll
import matplotlib.pyplot as plt
import numpy as np
import librosa as lb
import time

t0 = time.time()
folder = r'C:\Users\Loes\Documents\GitHub\Projectstage\audiovoorbeelden\heartbeat-sounds\set_a'
x, sr, duration = ll.retrieve_audio_data(folder, plot = False)
segmented_data = ll.segment_audio(x, duration)



xft = ll.compute_fft(segmented_data)

#  Extract features frequency domain
f = 1;
t = 1;
magnitude = []
phase = []

t1 = time.time()

for i in range(len(xft)):
    for f in range(len(xft[i])):
        for t in range(len(xft[i][f])):
            magnitude.append(np.abs(xft[i][f][t]))
            phase.append(np.angle(xft[i][f][t]))
            #print('magnitude = {} \nphase = {}'.format(magnitude, phase))

t2 = time.time()
print('first part: {:.3f}s \nsecond part: {:.3f}s \noverall: {:.3f}s'.format(t1-t0, t2-t1, t2-t0))



#a, X_list = ll.preprocess_data(folder, plot=False)

# =============================================================================
# n_training_samples = data.shape[0]
# n_dim = data.shape[1]
# 
# plt.figure()
# plt.plot(data[:,0], data[:,1], "bx")
# plt.show()
# =============================================================================

# Data plotten lukt niet, want de data heeft meer dan 2 dimenties.
# Kijk naar Xft en Xft_list
# Onset_detect doet het niet met Xft.
# Daarna verder anomaly detection uitvoeren mbv http://aqibsaeed.github.io/2016-07-17-anomaly-detection/