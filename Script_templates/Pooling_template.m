% Make sure that this function is set as you working directory, because the script will add the necessary dependencies.
clearvars
close all
 
gpuDevice; 
affine2d;
 
addpath(genpath('./helperfunctions'))
addpath(genpath('./helperfunctions/ext'))
 
% If the user didn’t save the processed data in mat format but has the scripts containing all the parameters, he can automatically run all the scripts corresponding to datasets to pool here.
 
fileNamesSciptToPool{1} = 'script1.m'; % Replace script1 by the name of your first file
fileNamesSciptToPool{2} = 'script2.m'; % Replace script2 by the name of your second file

% If you have more scripts corresponding to data to pool, add them here.
 
for i = 1:length(fileNamesSciptToPool)
     evalScript(fileNamesSciptToPool{i});
end
 
%%
% Now the pipeline will pool the replicates. 
clearvars
close all
 
fileNameBase = 'X'; % Enter instead of X the name common to all the files to pool (the user might have called the files for example: experiment_replicate1_Target.mat, experiment_replicate2_Target.mat, experiment_replicate1_DarkROIs.mat, experiment_replicate2_DarkROIs.mat. In this case fileNameBase = 'experiment')

% Now all the files with Target.mat at the end of their names will be pooled as data for target locations, and all the files with DarkROIs.mat at the end of their names will be pooled as data for dark locations. If the user didn’t use the same procedure in naming files, he should change the two following lines accordingly.
 
filesDark = dir([fileNameBase '*DarkROIs.mat']); 
filesExp = dir([fileNameBase '*Target.mat']);

for i = 1:length(filesExp)
    experimentFiles{i} = [filesExp(1).name];
    experimentFilesDarkROIs{i} = [filesDark(1).name];
end
 
h2 = createParentFigure;
htabgroup2 = uitabgroup(h2);
gepObj2 = guiExperimentPooling(htabgroup2,experimentFiles);
gepdrObj2 = guiExperimentPoolingDarkROI(htabgroup2,experimentFilesDarkROIs);
gdcorObj = guiDarkCorrector(htabgroup2,gepObj2,gepdrObj2);

% At this point, the user can look at fluorescence intensity time traces, plot rastergrams and perform correction for non-specific binding via the interface of the pipeline. 
