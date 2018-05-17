%Splitsen dataset in train en test set
clear, close, clc

workingdir = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\wavFiles\Dataset 1\';
base_folder = strcat(workingdir,'Normal\');
considerlogs = textread(strcat(workingdir,'filenames_wav.txt'), '%s');
mkdir(base_folder)
perc_train = 0.8;
perc_test = 0.2;

for i = 1:numel(considerlogs)
    selected_log = strcat(base_folder,considerlogs{i});
    fprintf('selected recording log: %s \n', selected_log)
    
    [x,Fs] = audioread(selected_log); %The input values from audioread() are dimensionless, scaled to -1<=x<1
    duration = floor(length(x)/Fs); %in seconds, normale dataset die nodig is, duurt 600 seconden.
    train_samples = perc_train*duration*Fs; %een aantal samples van x nemen voor trainen
    test_samples = perc_test*duration*Fs; %een aantal samples van x nemen voor testen
    train_set = x(1:train_samples);
    test_set = x(train_samples+1:train_samples+test_samples);  
    
    s1 = considerlogs{i};
    s2 = s1(1:end-8);
    %write test set
    write_name_test=strcat(s2,'_TEST','.wav'); % need to change for every file
    fullFileName_test = fullfile(base_folder, write_name_test);
    audiowrite(fullFileName_test,test_set,Fs,'BitsPerSample',16); % outputfile is normalized in apmlitude
    %write train set
    write_name_train=strcat(s2,'_TRAIN','.wav'); % need to change for every file
    fullFileName_train = fullfile(base_folder, write_name_train);
    audiowrite(fullFileName_train,train_set,Fs,'BitsPerSample',16); % outputfile is normalized in apmlitude
end

