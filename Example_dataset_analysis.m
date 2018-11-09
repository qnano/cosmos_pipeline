%%
close all
clearvars

addpath(genpath('.\helperfunctions'))

path = '.\Example_dataset\';

fileStructure.Grid =  [ path 'Grid.tif'];
fileStructure.Dark = [ path 'Dark.tif'];
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

gsdObj.paramsPreFilterFits.circularityMin=0.5;

gsdObj.paramsPreFilterFits.circularityMax=2;

gsdObj.paramsPreFilterFits.pH1Min=0;
gsdObj.paramsPreFilterFits.pH1Max=1;
gsdObj.paramsPreFilterFits.minPixelDist=7;
gsdObj.paramsPreFilterFits.clusterSizeMin=0;
gsdObj.paramsPreFilterFits.clusterSizeMax=100;

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

gadObj.spotsIncluded([2,5,6,10,12,15,16,19,32,35,44,45,55,64,69,74,77,...
    79,84,85,86,88,92,104,125,129,135,140,143,147,156,159,161,166,184,...
    187,194,201,203,214,238,242,265,267,269,279,280,284,288,297,298,301,...
    312,315,316,318,320,321,323,327,340,345,351,354,360,364,368,370,395,...
    396,398,414,418,423,426,429,432,442,443,445,446]) = 0;
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


%gdrObj.savevars([mfilename date 'DarkROIs.mat']);
%gadObj.savevars([mfilename date 'Target.mat']);

gadObj.prior.alpha = 1; 
gadObj.prior.kappa = 100; 
gadObj.prior.m = '[500 1000 1500 2000; 500 1000 1500 2000]'; 
gadObj.prior.v = 10; 
gadObj.prior.W = 10; 
gadObj.NTimes = 10; 
gadObj.NOrder= 4; 
gadObj.maximumlikelihood=false;
gadObj.bayesSignalFitGMM;


gadObj.prior.alpha = 1; 
gadObj.prior.kappa = 40; 
gadObj.prior.m = '[500 1000 1500 2000; 500 1000 1500 2000]'; 
gadObj.prior.v = 10; 
gadObj.prior.W = 10; 
gadObj.NTimes = 2; 
gadObj.NOrder= 4; 
gadObj.maximumlikelihood=true;
gadObj.bayesSignalFitGMM;

% % 
% % N=10;
% % [meanDwellTimeBS,meanTime2FirstEventBS] = gadObj.bootStrapRastergram(N); %takes N times randomly 90% of the data
% % 
% % sortState=1;
% % [meanDwellTimeBS,meanTime2FirstEventBS] = gadObj.bootStrapRastergram(N,sortState); %sort by state 1


% gdcorObj.getCorrectedOnrate('cdf');
% gdcorObj.getCorrectedOffrate('cdf');

