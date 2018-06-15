%Check normality of data
%% Histogram
%Show plot in window and also save it to a specified folder
close all
for ii= 1:numFeatures
    h = figure
    histogram(multi_feature_vectors_normaal(ii,:));
    xlabel(['Data of Feature ', num2str(ii)])
    ylabel('Frequency')
    title('Histogram')
    fname = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\images';
    filename= sprintf('Histogram feature %d',ii);
    saveas(h,fullfile(fname,filename), 'jpeg')
end
%% Normal plot
%Show plot in window and also save it to a specified folder
close all
for jj=1:numFeatures
    h = figure
    normplot(multi_feature_vectors_normaal(jj,:))
    xlabel(['Data of Feature ', num2str(jj)])
    fname = 'C:\Users\Gebruiker\Documents\GitHub\Projectstage\images';
    filename= sprintf('Normal plot feature %d',jj);
    saveas(h,fullfile(fname,filename), 'jpeg')
end    