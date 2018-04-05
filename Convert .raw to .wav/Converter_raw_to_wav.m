% Multiple file raw audio file conversion

list = dir ('test_raw_file.raw');

endianness='l';
fs=16000;
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