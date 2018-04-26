clear, close, clc
addpath('C:\Users\Loes\Documents\GitHub\Projectstage\Create Features')

reFs = 48000;
duration_window = 200;

path = 'C:\Users\Loes\Documents\GitHub\Projectstage\wavFiles\Waterkoker\';
filename_part_1 =  'mic_44100_s16le_channel_';
filename_part_2 = '_WAV.wav';
baseline_cellarray = BaselineLib(reFs, duration_window, path, filename_part_1, filename_part_2);