%%
close all
clearvars

addpath(genpath('.\helperfunctions'))

path = '.\Example_dataset\';

fileStructure.Grid = [path 'Grid.tif'];
fileStructure.Dark = [path 'Dark.tif'];
fileStructure.BeadData = [path 'Beads.tif'];
fileStructure.Data = [path 'Experiment.tif'];
fileStructure.rg = -1;
fileStructure.NFramesGroup = 1;
fileStructure.lambdal = 550;
fileStructure.lambdar = 647;
fileStructure.intT = 5;
fileStructure.flipCams = 0;

numberOfSummedFrames = 1;

h = createParentFigure;
cameraIndexDriftEst=1;

htabgroup = uitabgroup(h);
gdlObj = guiDataLoader(htabgroup,fileStructure,[],numberOfSummedFrames,1);
set(0, 'CurrentFigure', h)
gaObj = guiAlignmentw(htabgroup,gdlObj);
gdcObj = guiDriftCorr(htabgroup,gdlObj,cameraIndexDriftEst);
gsdObj = guiSpotDetector(htabgroup,gdlObj,[],[],[],[],true);
gadObj = guiHMMAnalyse(htabgroup,gdlObj,gaObj,gdcObj,gsdObj);
gdrObj = guiDarkROIsw(htabgroup,gdlObj,gaObj,gdcObj,gsdObj);
gdcorObj = guiDarkCorrector(htabgroup,gadObj,gdrObj);
gocObj = guiDataCollector(gdlObj,gaObj,gdcObj,gsdObj,gadObj,gdrObj);

gaObj.alignmentEst;

%%

C = [...
  235   299;...
  165    90];

gdcObj.driftEst(C);

% set paramsPreFilterFits
gsdObj.paramsPreFilterFits.circularityMin=0.8;
gsdObj.paramsPreFilterFits.circularityMax=1.4;
gsdObj.paramsPreFilterFits.pH1Min=0;
gsdObj.paramsPreFilterFits.pH1Max=1;
gsdObj.paramsPreFilterFits.minPixelDist=7;
gsdObj.paramsPreFilterFits.clusterSizeMin=0;
gsdObj.paramsPreFilterFits.clusterSizeMax=50;

% set paramsFilterFits 
gsdObj.paramsFilterFits.MinPhotons=50;
gsdObj.paramsFilterFits.MaxPhotons=Inf;
gsdObj.paramsFilterFits.MinBg=0;
gsdObj.paramsFilterFits.MaxBg=Inf;
gsdObj.paramsFilterFits.MinPValue=0;
gsdObj.paramsFilterFits.MaxPValue=1;
gsdObj.paramsFilterFits.MinPixelDist=2;
gsdObj.paramsFilterFits.MinCRLBSTD=0;
gsdObj.paramsFilterFits.MaxCRLBSTD=0.5;

gsdObj.detectSpots

gadObj.analyze

gadObj.spotsIncluded([]) = 0;
gadObj.rastergramStartframe = 80;

gadObj.setThreshold([150 1 1 4.6]);
gabDurVec = gadObj.getGabDurSliders;
gadObj.setGabDurSliders([5 5]);

gadObj.sortType = 3;
gadObj.sortState = 1;
gadObj.align = 0;
gadObj.delete = 0;
gadObj.rastergramStartframe = 80;
gadObj.colorArray = [1 1 1;0 0 1;1 0 0;1 0 0;0 1 0;1 0.4 0.4;0 0 0];
gadObj.numberOfSpotsForRastergram = [];
gadObj.runRastergram

gdrObj.analyze

gdrObj.spotsIncluded([]) = 0;

gabDurVec = gdrObj.getGabDurSliders;
gdrObj.setGabDurSliders([5 5]);

gdrObj.sortType = 3;
gdrObj.sortState = 1;
gdrObj.align = 0;
gdrObj.delete = 0;
gdrObj.rastergramStartframe = 80;
gdrObj.colorArray = [1 1 1;0 0 1;1 0 0;1 0 0;0 1 0;1 0.4 0.4;0 0 0];
gdrObj.numberOfSpotsForRastergram = [];
gdrObj.runRastergram

gdcorObj.getCorrectedOnrate('cdf');
gdcorObj.getCorrectedOffrate('cdf');

