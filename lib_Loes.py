import sys
import os
import numpy, scipy, matplotlib.pyplot as plt
import librosa, librosa.display


def retrieve_audio_data(folder, plot=True):
    files = []
    x_list = []
    for _, _, files in os.walk(folder):
<<<<<<< HEAD
        for filename in files[:5]:
=======
        for filename in files[:1]:
>>>>>>> c6584eec9bd07e77299e3e4c991df77818756772
            files.append(filename)
    #files = [filename for filename in files[:5] for _, _, files in os.walk(folder)]
    #Load the audio file simple_loop.wav into an array. With amplitude as values
    x, sr = librosa.load(os.path.join(folder, files[1]), sr=None)
<<<<<<< HEAD
    for filename in files[:5]:
=======
    for filename in files[:1]:
>>>>>>> c6584eec9bd07e77299e3e4c991df77818756772
        x, sr = librosa.load(os.path.join(folder, filename), sr=None) #sr is sampling rate
        x_list.append(x)
        duration = len(x) / float(sr)
        if plot == True:
            plt.figure()
            librosa.display.waveplot(x, sr=sr)
    return x_list, sr, duration

def segment_audio(data, duration, milliseconds = 200):
    ms = milliseconds
    n_segments = duration*1000 / float(ms)
    all_segments = []
    
    # how many frames in 1 millisecond? * ms
    frames_per_ms = int(len(data[0]) / (duration * 1000)) * ms
    for i in range(len(data)):
        segments = []
        from_frames_per_ms = 0
        to_frames_per_ms = 0
        for n in range(int(n_segments)):
            # begin at x[t=0] to x[t=time_segment]
            segments.append((data[0][from_frames_per_ms:to_frames_per_ms]))
            from_frames_per_ms = frames_per_ms
            to_frames_per_ms += frames_per_ms
        all_segments.append((segments))
    return all_segments

def compute_fft(segmented_data):
    xft_list = []
    xft_all = []
    for i in range(len(segmented_data)):
        for j in range(2,len(segmented_data[i])):
            xft = librosa.stft(segmented_data[i][j])
            xft_list.append(xft)
        xft_all.append(xft_list)
    return xft_all



# =============================================================================
# def segment_audio(x_all, sr):
#     frame_sz = int(0.100*sr)
#     segments_list = []
#     
#     for x in x_all:
#         #Find the times, in seconds, when onsets occur in the audio signal.
#         onset_frames = librosa.onset.onset_detect(x, sr=sr)
#         #Convert the onset frames into sample indices.
#         onset_samples = librosa.frames_to_samples(onset_frames)
#         
#         #----Segment the Audio----------
#         segments_list.append(numpy.array([x[i:i+frame_sz] for i in onset_samples]))
#     return segments_list
# 
# def concatenate_segments(segments_list, sr=22050, pad_time=0.300):
#     padded_segments = []
#     conc_signal = []
#     for segments in segments_list:
#         padded_segments.append([numpy.concatenate([segment, numpy.zeros(int(pad_time*sr))]) for segment in segments])
#     for padded_segment in padded_segments:
#         conc_signal.append(numpy.concatenate(padded_segment))
#     return conc_signal
# 
# def extract_features(segments_list, sr):
#     zero_crossing_rate_list = []
#     ind_list = []
#     segments_zcrs = []
#     for segments in segments_list:
#         #For each segment, compute the zero crossing rate.
#         zero_crossing_rate_list.append([sum(librosa.core.zero_crossings(segment)) for segment in segments])
#         
#     for zcrs in zero_crossing_rate_list:  
#         #Use argsort to find an index array, ind, such that segments[ind] is sorted by zero crossing rate.
#         ind_list.append(numpy.argsort(zcrs))
#         
#     for i in range(len(segments_list)):
#         segments_zcrs.append(segments_list[i][ind_list[i]])
#     
#     #Sort the segments by zero crossing rate, and concatenate the sorted segments.
#     concatenated_signal = concatenate_segments(segments_zcrs, sr)
#     return concatenated_signal
# 
# def preprocess_data(folder, plot=True):
#     x, sr, Xft_list = retrieve_audio_data(folder, plot=False)
#     segments_list = segment_audio(x, sr)
#     concatenated_signal = extract_features(segments_list, sr)
#     return concatenated_signal, Xft_list
# =============================================================================

