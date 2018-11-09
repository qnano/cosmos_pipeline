% Make sure that this function is set as you working directory, because the script will add the necessary dependencies.
close all
clearvars

% Add to path the folder containing required functions 
addpath(genpath('.\helperfunctions'))

% Here the user indicates where are located the dataset to be analyzed and how the experiment was done (number of frames per second, wavelengths of the two cameras etc). 
% Check gdlObj.fileStructure to get all single variables. For example, type in the command line “fileStructure.Grid” (without quotation marks) and copy-paste the answer.

fileStructure.Grid = ''; % write the answer between ''
fileStructure.Dark = ''; % write the answer between ''
fileStructure.BeadData = ''; % write the answer between ''
fileStructure.Data = ''; % write the answer between ''
fileStructure.rg = ; % write the answer before ; mark
fileStructure.NFramesGroup = ; % write the answer before ; mark
fileStructure.StartFrameSummation= ; % write the answer before ; mark
fileStructure.lambdal = ; % write the answer before ; mark
fileStructure.lambdar = ; % write the answer before ; mark
fileStructure.intT = ; % write the answer before ; mark
fileStructure.flipCams = ; % write the answer before ; mark

% Now the pipeline will load the datasets for the analysis. The user does not need to change anything in the code. 

h = createParentFigure;
cameraIndexDriftEst=1;

htabgroup = uitabgroup(h);
gdlObj = guiDataLoader(htabgroup,fileStructure);
set(0, 'CurrentFigure', h)
gaObj = guiAlignmentw(htabgroup,gdlObj);
gdcObj = guiDriftCorr(htabgroup,gdlObj,cameraIndexDriftEst);
gsdObj = guiSpotDetector(htabgroup,gdlObj,[],[],[],[],true);
gadObj = guiHMMAnalyse(htabgroup,gdlObj,gaObj,gdcObj,gsdObj);
gdrObj = guiDarkROIsw(htabgroup,gdlObj,gaObj,gdcObj,gsdObj);
gdcorObj = guiDarkCorrector(htabgroup,gadObj,gdrObj);
gocObj = guiDataCollector(gdlObj,gaObj,gdcObj,gsdObj,gadObj,gdrObj);

% Now the pipeline will perform alignment. 

gaObj.spotDetectorCam1.paramsFilterFits.minPixelDist = ; % write the answer before ; mark
gaObj.spotDetectorCam2.paramsFilterFits.minPixelDist = ; % write the answer before ; mark
gaObj.spotDetectorCam1.paramsPreFilterFits.clusterSizeMax = ; % write the answer before ; mark
gaObj.spotDetectorCam2.paramsPreFilterFits.clusterSizeMax = ; % write the answer before ; mark
gaObj.alignmentEst;

%%
 
% Now the pipeline will perform drift correction. 
% Type gdcObj.C to know which region was used for drift correction in you ongoing analysis and set region

C = [...
;... % Copy the first row of the answer here before ; mark
]; % Copy the second row of the answer here before ; mark
 
gdcObj.driftEst(C);
 
% Now the pipeline will detect locations (called here spits) of target molecules and mobile components. 

% set parameters used for pre-filtering spots

gsdObj.paramsPreFilterFits.circularityMin = ; % write the answer before ; mark
gsdObj.paramsPreFilterFits.circularityMax = ; % write the answer before ; mark
gsdObj.paramsPreFilterFits.PH1Min = ; % write the answer before ; mark
gsdObj.paramsPreFilterFits.PH1Max = ; % write the answer before ; mark
gsdObj.paramsPreFilterFits.minPixelDist = ; % write the answer before ; mark
gsdObj.paramsPreFilterFits.clusterSizeMin = ; % write the answer before ; mark
gsdObj.paramsPreFilterFits.clusterSizeMax = ; % write the answer before ; mark
 
% Set parameters used for filtering spots

gsdObj.paramsFilterFits.MinPhotons = ; % write the answer before ; mark
gsdObj.paramsFilterFits.MaxPhotons = ; % write the answer before ; mark
gsdObj.paramsFilterFits.MinBg = ; % write the answer before ; mark
gsdObj.paramsFilterFits.MaxBg = ; % write the answer before ; mark
gsdObj.paramsFilterFits.MinPValue = ; % write the answer before ; mark
gsdObj.paramsFilterFits.MaxPValue = ; % write the answer before ; mark
gsdObj.paramsFilterFits.MinPixelDist = ; % write the answer before ; mark
gsdObj.paramsFilterFits.MinCRLBSTD = ; % write the answer before ; mark
gsdObj.paramsFilterFits.MaxCRLBSTD = ; % write the answer before ; mark
 
gsdObj.detectSpots;

% Now the pipeline will plot the fluorescence intensity time traces of target molecules and mobile components and detect co-localization events. 

% Set the criteria for co-localization. First, type “gadObj.getThreshold” (without quotation marks) in the command line to get the thresholds used in the ongoing analysis. Second, type “gadObj.getGabDurSliders” (without quotation marks) in the command line to get the values of gap closing and minimum value of a co-localization event.
gadObj.setThreshold([X X X X]); % replace Xs by the values from gadObj.getThreshold
gabDurVec = gadObj.getGabDurSliders;
gadObj.setGabDurSliders([X X]); % replace Xs by the values from gadObj.getGabDurSliders


gadObj.analyze
 
% To find indices of spots, which the user decided to exclude from the analysis, type in the command line “a=find(~gadObj.spotsIncluded)” (without the quotation marks). Then transpose the vertical vector to horizontal vector by typing “b=a.'” (without the quotation marks).
gadObj.spotsIncluded([]) = 0; % Enter the values from b between []. If no spot was excluded from the analysis, leave this line without modification.

% Now the pipeline will plot the rastergram summarizing the co-localization events. 

% Enter rastergram settings
gadObj.rastergramStartframe = ; % write the answer before ; mark
gadObj.sortType = ; %none=0;duration of first=1;duration of last=2;binding of first=3;
gadObj.sortState = 1; % state to sort detult is 1
gadObj.align = ; %none=0;first binding event=1;last departure event=2;
gadObj.delete = ; %true = 1 and false = 0
gadObj.rastergramStartframe = ; %frame to start rastergram
gadObj.colorArray = ; %rgb values of the color of each state e.g [1 1 1;0 0 1;1 0 0;1 0 0;0 1 0;1 0.4 0.4;0 0 0]
gadObj.numberOfSpotsForRastergram = []; % empty for using all states


gadObj.runRastergram
 
% Now the pipeline will analyze non-specific binding at dark locations. All the parameters (the criteria for co-localization and the values of gap closing and minimum value of a co-localization event) are automatically set at the same values than for analysis of binding at the target molecules.
 
gdrObj.analyze
 
% To find indices of spots, which the user decided to exclude from the analysis, type in the command line “c=find(~gdrObj.spotsIncluded)” (without the quotation marks). Then transpose the vertical vector to horizontal vector by typing “d=c.'” (without the quotation marks).

gdrObj.spotsIncluded([]) = 0; % Enter the values from d between []. If no spot was excluded from the analysis, leave this line without modification.

gabDurVec = gdrObj.getGabDurSliders;
gdrObj.setGabDurSliders([5 5]);

% Now the pipeline will plot the rastergram summarizing the co-localization events. All rastergram settings are automatically set at the same values than for analysis of binding at the target molecules.

% Copy rastergarm settings from the analysis module (gadObj)
gdrObj.copyRastergramSettings(gadObj);

gdrObj.rastergramStartframe = 80;

gdrObj.runRastergram

% Now the pipeline saves all processed data (alignment of cameras, drift correction, fluorescence intensity time traces, and segmentation into binding events) as .mat files. This is done for target molecules and dark locations separatly.

gdrObj.savevars([mfilename date 'DarkROIs.mat']);
gadObj.savevars([mfilename date 'Target.mat']);

% Set HMM settings
gadObj.prior.alpha = ; Prior \alpha
gadObj.prior.kappa = ; Prior \kappa
gadObj.prior.m = ; Prior m
gadObj.prior.v = ; Prior v
obj.prior.W = ; Prior W^{1/2}
gadObj.NTimes = ; N Samples
gadObj.NOrder= ; Maximum Order
