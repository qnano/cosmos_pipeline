clearvars
close all

gpuDevice; 
affine2d;

addpath(genpath('./helperfunctions'))
addpath(genpath('./helperfunctions/ext'))


%% Pooling example
clearvars
close all

path = './Example_processed_data/';

fileNameBase = [path 'Example_rep']

filesDark = dir([fileNameBase '*DarkROIs.mat']);
filesExp = dir([fileNameBase '*Target.mat']);
for i = 1:length(filesExp)
    experimentFiles{i} = [path filesExp(i).name];
    experimentFilesDarkROIs{i} = [path filesDark(i).name];
end

h2 = createParentFigure;
htabgroup2 = uitabgroup(h2);
gepObj2 = guiExperimentPooling(htabgroup2,experimentFiles);
gepdrObj2 = guiExperimentPoolingDarkROI(htabgroup2,experimentFilesDarkROIs);
gdcorObj = guiDarkCorrector(htabgroup2,gepObj2,gepdrObj2);
