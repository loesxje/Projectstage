% Multiple file raw audio file conversion
clear, close, clc
path = 'C:\Users\Loes\Documents\GitHub\Projectstage\Convert .raw to .wav\';
k = 7;
for i = 0:k
    list(i) = sprintf('mic_44100_s16le_channel_%d_.raw',i)
    [x,Fs] = audioread(strcat(path,variable));
end



%list = dir ('mic_44100_s16le_channel_8.raw');

endianness='l';
fs=44100;
datatype='int16';
for i = 1 : length(list)
    filename=list(i).name;
    fid_R = fopen(filename,'r',endianness);
    data = fread(fid_R,'int16'); % Raw data of input audio file. 
    fclose(fid_R); 
% %     o=audioplayer(data,16000);
% %     play(o)
    data=int16(data);
    s1=list(i).name;
    s2=s1(1:end-4);
    write_name=strcat(s2,'_WAV.wav'); % need to change for every file
    audiowrite(write_name,data,fs,'BitsPerSample',16); % outputfile is normalized in apmlitude   
end