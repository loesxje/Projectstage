clear, close, clc

workingdir = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\';
basefolderlogs = strcat(workingdir,'rawFiles\');
target_folder = strcat(workingdir,'wavFiles\');
considerlogs = textread(strcat(workingdir,'rawFiles\filenames.txt'), '%s');
mkdir(target_folder)

endianness='l';
fs=44100;
datatype='int16';

for i = 1:numel(considerlogs)
    selected_log = strcat(basefolderlogs,considerlogs{i});
    fprintf('selected recording log: %s \n', selected_log)
    fid_R = fopen(selected_log,'r',endianness);
    data = fread(fid_R,'int16'); % Raw data of input audio file. 
    fclose(fid_R); 
    data=int16(data);
    s1 = considerlogs{i};
    s2 = s1(1:end-4);
    write_name=strcat(s2,'_WAV.wav'); % need to change for every file
    fullFileName = fullfile(target_folder, write_name);
    audiowrite(fullFileName,data,fs,'BitsPerSample',16); % outputfile is normalized in apmlitude
end