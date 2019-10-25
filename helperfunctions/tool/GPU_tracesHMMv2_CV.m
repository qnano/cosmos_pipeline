close all
clearvars

%%
addpath(genpath('./helperfunctions'))
addpath(genpath('./rasterscripts'))

%% Load alignment data
% From file
% fileStructure(1).Grid = '/2016_0203_CamCal/Grid_100_120_1fps_2.tif';
% fileStructure(1).Dark = '/2016_0203_CamCal/Dark_100_120_1fps.tif';
% fileStructure(1).BeadData = '/2016_0503_Cas9/beads2_map.tif';
% fileStructure(1).Data = '/2016_0503_Cas9/merged004.tif';
% fileStructure(1).Grid = '/2016_0203_CamCal/Grid_100_120_1fps_2.tif';
% fileStructure(1).Dark = '/2016_0203_CamCal/Dark_100_120_1fps.tif';
% fileStructure(1).BeadData = '/2016_0406_Cas9/beads2_map.tif';
% fileStructure(1).Data = '/2016_0406_Cas9/merged002.tif';
% fileStructure(1).BeadData = '/2016_0203_TtAgo/beads2_map.tif';
% fileStructure(1).Data = '/2016_0203_TtAgo/merged007.tif';


fileStructure(1).Grid = './Karina_20160511/Grid_90_90_5fps_1.tif';
fileStructure(1).Dark = './Karina_20160511/Dark_90_90_5fps.tif';
fileStructure(1).BeadData = './Wes_Victor/Beads_1268/beads2_map.tif';
fileStructure(1).Data = './Wes_Victor/1268_merged002_perfect.tif';
%% Callibration Data
a = bfopen( [pwd  fileStructure(1).Grid]);
geller = double(cell2mat(permute(a{1}(1:min(100,end),1),[3 2 1])));

a = bfopen( [pwd  fileStructure(1).Dark] );
bg = double(cell2mat(permute(a{1}(1:min(100,end),1),[3 2 1])));

gellerCam2 = geller(1:obj.fileStructure.imageDimensions(1),...
    1:obj.fileStructure.imageDimensions(2),:); % swapped Cam1 and Cam2 - VS
gellerCam1 = geller(1:obj.fileStructure.imageDimensions(1),...
    obj.fileStructure.imageDimensions(2)+1:2*obj.fileStructure.imageDimensions(2),:);

bgCam2 = bg(1:obj.fileStructure.imageDimensions(1),1:obj.fileStructure.imageDimensions(2),:); % swapped Cam1 and Cam2 - VS
bgCam1 = bg(1:obj.fileStructure.imageDimensions(1),...
    obj.fileStructure.imageDimensions(2)+1:2*obj.fileStructure.imageDimensions(2),:);

outCam1 = cal_readnoise(gellerCam1, bgCam1);
outCam2 = cal_readnoise(gellerCam2, bgCam2);

PSFSigma=1.39;
nFrames = 10;

a = bfopen( [pwd  fileStructure(1).BeadData] );
data = double(cell2mat(permute(a{1}(1:min(end,nFrames),1),[3 2 1])));
data=(data-repmat(mean(bg,3),[1 1 size(data,3)]));

beadDataCam2 = data(1:obj.fileStructure.imageDimensions(1),...
    1:obj.fileStructure.imageDimensions(2),1:min(end,nFrames))*outCam1(2); % swapped Cam1 and Cam2 - VS
beadDataCam1 = data(1:obj.fileStructure.imageDimensions(1),...
    obj.fileStructure.imageDimensions(2)+1:2*obj.fileStructure.imageDimensions(2),...
    1:min(end,nFrames))*outCam2(2);

%%
PSFSigma=1.39;
nFrames = 1e4;

a = bfopen( [pwd  fileStructure(1).Data] );
data = double(cell2mat(permute(a{1}(1:min(end,nFrames),1),[3 2 1])));
data=(data-repmat(mean(bg,3),[1 1 size(data,3)]));

dataCam2 = data(1:obj.fileStructure.imageDimensions(1),...
    1:obj.fileStructure.imageDimensions(2),1:min(end,nFrames))*outCam1(2); % swapped Cam1 and Cam2 - VS
dataCam1 = data(1:obj.fileStructure.imageDimensions(1),...
    obj.fileStructure.imageDimensions(2)+1:2*obj.fileStructure.imageDimensions(2),...
    1:min(end,nFrames))*outCam2(2);

%% GLRT Based Detection

[coordsCam1,detParCam1,cutProcessCam1] = LLRMapv2(beadDataCam1,PSFSigma);
[coordsCam2,detParCam2,cutProcessCam2] = LLRMapv2(beadDataCam2,PSFSigma);

%% Pre Filter Detection Clusters Cam1
paramsPreFilterFits = getDefaultParamsPreFilterFits;
paramsPreFilterFits.minPixelDist=2;
paramsPreFilterFits.clusterSizeMax=25;
paramsPreFilterFits.circularityMax=2;
[ maskPreFiltCam1 ] =  preFilterFits(coordsCam1,detParCam1,paramsPreFilterFits);

% optionsLLR = dipSubLoc2DSetOptions;
% optionsLLR.BoxCenters = [coordsCam1(find(maskPreFiltCam1),1) coordsCam1(find(maskPreFiltCam1),2) coordsCam1(find(maskPreFiltCam1),3)];
% optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
% C=jet(256); 
% index = double(floor(255.*(detParCam1.PH1(:,1)))+1);
% optionsLLR.BoxColor = C(index,:);
% optionsLLR.plotBoxes = 1;
% optionsLLR.im = dip_image(stretch(permute( single(squeeze(stretch(cut(beadDataCam1,[size(cutProcessCam1,2) size(cutProcessCam1,2) size(cutProcessCam1,3)])))),[1 2 3])),'uint8')';
% 
% h = dipSubLoc2D(optionsLLR);  

%% Pre Filter Detection Clusters Cam2
paramsPreFilterFits = getDefaultParamsPreFilterFits;
paramsPreFilterFits.minPixelDist=2;
paramsPreFilterFits.clusterSizeMax=25;
paramsPreFilterFits.circularityMax=2;
[ maskPreFiltCam2 ] =  preFilterFits(coordsCam2,detParCam2,paramsPreFilterFits);

% optionsLLR = dipSubLoc2DSetOptions;
% optionsLLR.BoxCenters = [coordsCam2(maskPreFiltCam2,1) coordsCam2(maskPreFiltCam2,2) coordsCam2(maskPreFiltCam2,3)];
% optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
% C=jet(256); 
% index = double(floor(255.*(detParCam2.PH1(maskPreFiltCam2,1)))+1);
% optionsLLR.BoxColor = C(index,:);
% optionsLLR.plotBoxes = 1;
% optionsLLR.im = dip_image(stretch(permute( single(squeeze(stretch(cut(beadDataCam2,[size(cutProcessCam2,2) size(cutProcessCam2,2) size(cutProcessCam2,3)])))),[1 2 3])),'uint8')';
% h = dipSubLoc2D(optionsLLR);  

%% MLE Fit Intensities Cam1

paramsFit = getDefaultParamsFit;
coodsUnCut=round(coordsCam1+(1.5*(2*PSFSigma+1)-0.5).*[ones(size(coordsCam1,1),1) ones(size(coordsCam1,1),1) zeros(size(coordsCam1,1),1)]);
paramsFit.FitSigma=true;

[ rawFitResultsCam1 ] = fitBoxCenters( single(squeeze(beadDataCam1)),coodsUnCut,paramsFit);

%%
paramsFilterFits = getDefaultParamsFilterFits;
paramsFilterFits.minPixelDist=2;


[ maskFilt1 ] =  filterFits(rawFitResultsCam1,paramsFilterFits);


% optionsLLR = dipSubLoc2DSetOptions;
% optionsLLR.BoxCenters = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1) rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,2) rawFitResultsCam1.Frame(maskPreFiltCam1&maskFilt1)];
% optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
% C=jet(256); 
% index = double(floor(255.*(detParCam1.PH1(maskPreFiltCam1&maskFilt1,1)))+1);
% optionsLLR.BoxColor = C(index,:);
% optionsLLR.plotBoxes = 1;
% optionsLLR.im = beadDataCam1; %dip_image(stretch(permute( single(squeeze(stretch(cut(dataCam1,[size(cutProcessCam1,2) size(cutProcessCam1,2) size(cutProcessCam1,3)])))),[1 2 3])),'uint8')';
% h = dipSubLoc2D(optionsLLR);  

%%

paramsFit = getDefaultParamsFit;
coodsUnCut=round(coordsCam2+(1.5*(2*PSFSigma+1)-0.5).*[ones(size(coordsCam2,1),1) ones(size(coordsCam2,1),1) zeros(size(coordsCam2,1),1)]);
paramsFit.FitSigma=true;

[ rawFitResultsCam2 ] = fitBoxCenters( single(squeeze(beadDataCam2)),coodsUnCut,paramsFit);

%%
paramsFilterFits = getDefaultParamsFilterFits;
paramsFilterFits.minPixelDist=2;
[ maskFilt2 ] =  filterFits(rawFitResultsCam2,paramsFilterFits);

% optionsLLR = dipSubLoc2DSetOptions;
% optionsLLR.BoxCenters = [rawFitResultsCam2.Coord(maskPreFiltCam2&maskFilt2,1) rawFitResultsCam2.Coord(maskPreFiltCam2&maskFilt2,2) rawFitResultsCam2.Frame(maskPreFiltCam2&maskFilt2)];
% optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
% C=jet(256); 
% index = double(floor(255.*(detParCam2.PH1(maskPreFiltCam2&maskFilt2,1)))+1);
% optionsLLR.BoxColor = C(index,:);
% optionsLLR.plotBoxes = 1;
% optionsLLR.im = beadDataCam2; 
% h = dipSubLoc2D(optionsLLR);  

%%
mGridCam1 = mean(beadDataCam1,3);
mGridCam2 = mean(beadDataCam2,3);
sv2 = findshift(mGridCam1,mGridCam2,'iter');

dipshow(joinchannels('RGB',squeeze(shift(mGridCam1,-sv2)),squeeze(mGridCam2)))
[out,R] = find_affine_trans(mGridCam1, mGridCam2, [1, 1, sv2', 0])

[image_out,R] = affine_trans(mGridCam1, out(1:2), out(3:4), out(5));
dipshow(joinchannels('RGB',squeeze(mGridCam2),squeeze(image_out)))


% 
cpPoints = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1) rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,2) rawFitResultsCam1.Frame(maskPreFiltCam1&maskFilt1)];
cpPointsMoved  = [rawFitResultsCam2.Coord(maskPreFiltCam2&maskFilt2,1) rawFitResultsCam2.Coord(maskPreFiltCam2&maskFilt2,2) rawFitResultsCam2.Frame(maskPreFiltCam2&maskFilt2)];
tformNoised = affine2d(R');
cpPointsUnMoved = transformPointsInverse(tformNoised,cpPoints(:,1:2));

[idx, disttFormEst] = knnsearch(cpPointsMoved(:,1:2),cpPointsUnMoved);

im1 = smooth(coord2image(cpPointsMoved(idx(disttFormEst<1),1:2),...
    [obj.fileStructure.imageDimensions(1) obj.fileStructure.imageDimensions(2)]),2);
im2b = smooth(coord2image(cpPointsUnMoved(disttFormEst<1,:),...
    [obj.fileStructure.imageDimensions(1) obj.fileStructure.imageDimensions(2)]),2);
dipshow(joinchannels('RGB',squeeze(im1),squeeze(im2b)))


a = pinv([cpPointsMoved(idx(disttFormEst<1),1:2) ones(size(cpPointsMoved(idx(disttFormEst<1),:) ,1),1)])* [cpPointsUnMoved(disttFormEst<1,:) ones(size(cpPointsUnMoved(disttFormEst<1,:),1),1)];

Total = a*(R');
Total(:,3)=[0;0;1;];
tformTotal = affine2d(Total);
cpPointsUnMoved = transformPointsInverse(tformTotal,cpPoints(:,1:2));
 
[~, disttFormEst] = knnsearch(cpPointsMoved(idx(disttFormEst<1),1:2),transformPointsInverse(tformTotal,cpPoints(disttFormEst<1,1:2) ));
% close all

scatter3(cpPointsMoved(idx(disttFormEst<1),1),cpPointsMoved(idx(disttFormEst<1),2),cpPointsMoved(idx(disttFormEst<1),3),20,disttFormEst(disttFormEst<1),'filled')
colormap(jet);
caxis manual
caxis([0 max(1,max(disttFormEst))]);
colorbar;
shading interp;
% % smooth(coord2imagev2(cpPoints(idx,:),sensorSize,dist),10)
fprintf('Experimental mean square error: %0.2g\n',mean(disttFormEst.^2)) % estimated transform

%%
dataCam1_moved=dataCam1;

for ww=1:length(dataCam1(1,1,:))
    [dataCam1_moved(:,:,ww),R1] = affine_trans(dataCam1(:,:,ww), out(1:2), out(3:4), out(5));
end
% dataCam1_moved=dataCam1_moved*50;

a=joinchannels('RGB', stretch(dataCam1_moved(:,:,1:600),93,100,0,10000), stretch(dataCam2(:,:,1:600),85,100,0,7000));
dipshow(a)

%%
e=uint8(dip_array(a));
filename='2016_0503_Cas9_exp004';
VidObj=VideoWriter(filename,'Uncompressed AVI');
VidObj.FrameRate = 15;
% VidObj.Quality=100;
% VidObj.LosslessCompression=true;

open(VidObj);

for ii=1:size(e,3)
    curr_fr=squeeze(e(:,:,ii,:));
    curr_crop=curr_fr(250:500,50:463,:);
    writeVideo(VidObj,curr_crop);
end
close(VidObj);

%% plot pixel drift shift 
% C = [...
%    207   200;...
%     67    68];
% 
% if isempty(C)
%     h = dipshow(dataCam1(:,:,1));
%     [~,C] = dipcrop(h);
% % close(h)
% end

h = dipshow(dataCam1(:,:,1));
[~,C] = dipcrop(h);

cutRegion = double(dataCam1(C(1,1):C(1,1)+C(2,1),C(1,2):C(1,2)+C(2,2),:));

for i=1:size(dataCam1,3)-1
%     sv(:,i) = findshift(cutProcessCam2(ROI(1,:),ROI(2,:),i),cutProcessCam2(ROI(1,:),ROI(2,:),i+1),'iter');
    sv(:,i) = findshift(cutRegion(:,:,i),cutRegion(:,:,i+1),'iter');
end

figure
subplot(1,2,1)
plot(sv','linewidth',2)
legend('x-drift','y-drift')

ylabel('Drift [pixel]')
xlabel('t [frames]')
grid on

subplot(1,2,2)
plot(cumsum(sv'),'linewidth',2)
legend('x-drift','y-drift')

ylabel('Cumultive Drift [pixel]')
xlabel('t [frames]')
grid on

drift=cumsum(sv'); drift=[[0 0]; drift];
 
%% Selcect ROIs for dynamic estimation
% GLRT Based Detection
NFramesDetect = 2;

[canCoordsCam1,detParCam1,cutProcessCam1] = LLRMapv2(dataCam1(:,:,1:NFramesDetect),PSFSigma,0,0.5,0.05,8,true);

figure
subplot(1,3,1)
hist(detParCam1.circularity)
subplot(1,3,2)
hist(detParCam1.PH1)
subplot(1,3,3)
hist(detParCam1.clusterSize)

% Pre Filter Detection Clusters Cam1
paramsPreFilterFits = getDefaultParamsPreFilterFits;
paramsPreFilterFits.minPixelDist=3;
paramsPreFilterFits.clusterSizeMax=20;
paramsPreFilterFits.circularityMax=1.2;
[ maskPreFiltCam1 ] =  preFilterFits(canCoordsCam1,detParCam1,paramsPreFilterFits);

% optionsLLR = dipSubLoc2DSetOptions;
% optionsLLR.BoxCenters = [coordsCam1(find(maskPreFiltCam1),1) coordsCam1(find(maskPreFiltCam1),2) coordsCam1(find(maskPreFiltCam1),3)];
% optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
% C=jet(256); 
% index = double(floor(255.*(detParCam1.PH1(:,1)))+1);
% optionsLLR.BoxColor = C(index,:);
% optionsLLR.plotBoxes = 1;
% optionsLLR.im = dip_image(stretch(permute( single(squeeze(stretch(cut(dataCam1(:,:,1:NFramesDetect),[size(cutProcessCam1,2) size(cutProcessCam1,2) NFramesDetect])))),[1 2 3])),'uint8')';
% 
% h = dipSubLoc2D(optionsLLR);  

% MLE Fit Intensities Cam1

paramsFit = getDefaultParamsFit;
coodsUnCut=round(canCoordsCam1+(1.5*(2*PSFSigma+1)-0.5).*[ones(size(canCoordsCam1,1),1) ones(size(canCoordsCam1,1),1) zeros(size(canCoordsCam1,1),1)]);
paramsFit.FitSigma=true;

[ rawInitialFitResultsCam1 ] = fitBoxCenters( single(squeeze(dataCam1)),coodsUnCut,paramsFit);

%
paramsFilterFits = getDefaultParamsFilterFits;
paramsFilterFits.minPixelDist=7;


[ maskFilt1 ] =  filterFits(rawInitialFitResultsCam1,paramsFilterFits);

% 
% optionsLLR = dipSubLoc2DSetOptions;
% optionsLLR.BoxCenters = [rawInitialFitResultsCam1.Coord(maskPreFiltCam1 & maskFilt1,1) rawInitialFitResultsCam1.Coord(maskPreFiltCam1 & maskFilt1,2) rawInitialFitResultsCam1.Frame(maskPreFiltCam1 & maskFilt1)];
% optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
% C=jet(256); 
% index = double(floor(255.*(detParCam1.PH1(maskPreFiltCam1& maskFilt1,1)))+1);
% optionsLLR.BoxColor = C(index,:);
% optionsLLR.plotBoxes = 1;
% optionsLLR.im = dataCam1; %dip_image(stretch(permute( single(squeeze(stretch(cut(dataCam1,[size(cutProcessCam1,2) size(cutProcessCam1,2) size(cutProcessCam1,3)])))),[1 2 3])),'uint8')';
% h = dipSubLoc2D(optionsLLR);  

%%

coords = rawInitialFitResultsCam1.Coord(maskPreFiltCam1 & maskFilt1,:);
frames = rawInitialFitResultsCam1.Frame(maskPreFiltCam1 &maskFilt1,:);

photons1 = rawInitialFitResultsCam1.Photons(maskPreFiltCam1 & maskFilt1,:);
coords_1st=coords(frames==0,:);


traces1_allfr=zeros(size(coords_1st,1),NFramesDetect);
% traces2_allfr=zeros(size(coords_1st,1),NFramesDetect);

for w=1:NFramesDetect
    coords_nth=coords_1st-repmat(drift(w,:),size(coords_1st,1),1);
    coords_mapped=transformPointsInverse(tformTotal,coords_nth);
    [indx1,dist1]=knnsearch(coords(frames==w-1,:),coords_nth);
    
    % traces2=zeros(size(coords_corr,1),1);
    photons1_fr=photons1(frames==w-1);
    traces1=photons1_fr(indx1); traces1(dist1>1.5)=NaN;
    traces1_allfr(:,w)=traces1; % each row is a trace (Cam1)
end
% coords = rawInitialFitResultsCam1.Coord(maskPreFiltCam1 & maskFilt1,:);
% frames = rawInitialFitResultsCam1.Frame(maskPreFiltCam1 & maskFilt1,:);

%%
segmentedSpots = coords(~isnan(sum(traces1_allfr,2)),:);

% C =[20    100;...
%    500   500];
C=[];
if isempty(C)
    h=dipshow(joinchannels('RGB',squeeze(mGridCam2),squeeze(image_out)),'lin');
    [~,C] = dipcrop(h);

    close(h)
end

tSegmentedSpots = transformPointsInverse(tformTotal,segmentedSpots);
idx = C(1,1) < tSegmentedSpots(:,1) &...
C(1,1)+C(2,1) > tSegmentedSpots(:,1)&...
C(1,2) < tSegmentedSpots(:,2) &...
C(1,2)+C(2,2) > tSegmentedSpots(:,2);


dipshow(dataCam1(:,:,1))
hold on
plot(tSegmentedSpots(idx,1),tSegmentedSpots(idx,2),'xg')
hold on
plot(segmentedSpots(idx,1),segmentedSpots(idx,2),'xr')

segmentedSpots = segmentedSpots(idx,:);

%%
coordsCam1=[];
coordsCam2=[];
for w=1:size(dataCam1,3)
    coords_nth = segmentedSpots-repmat(drift(w,:),size(segmentedSpots,1),1);
    coordsCam1=cat(1,coordsCam1,[coords_nth (w-1).*ones(size(segmentedSpots,1),1)]);
    coordsCam2=cat(1,coordsCam2,[transformPointsInverse(tformTotal,coords_nth) (w-1).*ones(size(segmentedSpots,1),1)]); 
end

% optionsLLR = dipSubLoc2DSetOptions;
% optionsLLR.BoxCenters = coordsCam1;
% optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
% optionsLLR.BoxColor = [1 0 0]; 
% optionsLLR.plotBoxes = 1;
% optionsLLR.im = dataCam1; 
% h = dipSubLoc2D(optionsLLR);  
% 
% optionsLLR = dipSubLoc2DSetOptions;
% optionsLLR.BoxCenters = coordsCam2;
% optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
% optionsLLR.BoxColor = [1 0 0];
% optionsLLR.plotBoxes = 1;
% optionsLLR.im = dataCam2; %dip_image(stretch(permute( single(squeeze(stretch(cut(dataCam1,[size(cutProcessCam1,2) size(cutProcessCam1,2) size(cutProcessCam1,3)])))),[1 2 3])),'uint8')';
% h = dipSubLoc2D(optionsLLR);  

% MLE Fit Intensities Cam1
paramsFit = getDefaultParamsFit;
paramsFit.Iterations=16;
paramsFit.FitSigma=true;

[ rawFitResultsCam1 ] = fitBoxCenters( single(squeeze(dataCam1)),coordsCam1,paramsFit);


[ rawFitResultsCam2 ] = fitBoxCenters( single(squeeze(dataCam2)),coordsCam2,paramsFit);




%% Just show 1 spot in camera data2
% timeSteps = 20;
% 
% 
% h = dipshow(dataCam2);
% ha = findall(h,'type','axes');
% 
% diptruesize(h,200)
% movieFileName2 = 'spots.avi';
% 
% % Prepare VideoWriter object
% TimeStep=23;
% vidObj = VideoWriter(movieFileName2);
% vidObj.FrameRate = round(TimeStep);
% 
% %open avi file
% open(vidObj);
% Nspot = randi(size(photons1,1)/size(dataCam2,3),1,1);
% for ii = 0:timeSteps:size(dataCam2,3)-1
%     dipmapping(h,'slice',ii-1)
% %     h3= title(sprintf('spot ID = %d',Nspot))
%     hold on
%     idx = find(frames2==ii);
%     h1=plot(coordsCam2(idx(Nspot),1),coordsCam2(idx(Nspot),2),'xg');
%     h2= plot(coords2(idx(Nspot),1),coords2(idx(Nspot),2),'xb');
%     % Write each frame to the file.
%     currFrame = getframe(ha);
%     data = uint8(extend(currFrame.cdata,[194 194 3]));
%     currFrame.cdata = data;
%     writeVideo(vidObj,currFrame);
%     delete(h1);
%     delete(h2);
% %     delete(h3);
% end
% 
% %close avi file
% close(vidObj);
% close(h);
% 
% 
% figure
% 
% plot(photons2I(Nspot,:))
% hold on
% plot(bg2I(Nspot,:))
% plot(traces2_allfr(Nspot,:)*max(photons2I(Nspot,:)))
% legend('Intensity','Background','Binary')

%%
% This will make an Interval Data Structure from the traces (traces1_allfr and traces2_allfr) generated from 
% the output of camCal.m
% movieFileName2 = 'spotsHMM.avi';
% [traces2_allfr, seq, transmat1, prior1, mu1, Sigma1, mixmat1] =  HMMKEstimate(rawFitResultsCam1,rawFitResultsCam2,coordsCam2,tformTotal);
% 1/transmat1(1,2)
% 1/transmat1(2,1)



%%
close all
vid=[0:size(traces2_allfr,2)-1];

traces2_allfr_bin=traces2_allfr;
traces2_allfr_bin(traces2_allfr_bin>0)=1;
traces2_allfr_bin(:,end)=[]; % because time refers to the start of the frame, need to do this so that each frame has start time and end time 
traces2_allfr_bins=num2str(traces2_allfr_bin); 
traces2_allfr_binsc=num2cell(traces2_allfr_bins,2); % cell array is needed for regexp to work
traces2_allfr_binsc=regexprep(traces2_allfr_binsc,' ',''); % this makes a searchable cell array where each cell is a binary trace in string format



Intervals.CumulativeIntervalArray=[]; %initialize array

for i=1:length(traces2_allfr_binsc)
% Key (Larry's flags):
% "1..." = -3
% "...1" = 3
% "0..." = -2
% "...0" = 2
% "1...0...1" = 0
% "0...1...0" = 1
exp=traces2_allfr_binsc{i,1};
[start_flag_3, end_flag_3]=regexp(exp,'^1+');
[start_flag_2, end_flag_2]=regexp(exp,'^0+');
[start_flag0, end_flag0]=regexp(exp,'(?<=1)0{2,}(?=1)'); % at least 2 zeros flanked by ones, 1 frame events are discarded
[start_flag1, end_flag1]=regexp(exp,'(?<=0)1{2,}(?=0)'); % at least 2 ones flanked by zeros, 1 frame events are discarded
[start_flag2, end_flag2]=regexp(exp,'0+$');
[start_flag3, end_flag3]=regexp(exp,'1+$');
if ~isempty(start_flag_3)
    Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [-3 start_flag_3 end_flag_3 ...
        end_flag_3-start_flag_3+1 vid(end_flag_3+1)-vid(start_flag_3) 0 i]];
end
if ~isempty(start_flag_2)
    Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [-2 start_flag_2 end_flag_2 ...
        end_flag_2-start_flag_2+1 vid(end_flag_2+1)-vid(start_flag_2) 0 i]];
end
if ~isempty(start_flag0)
    Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [zeros(length(start_flag0),1) ...
        start_flag0' end_flag0' end_flag0'-start_flag0'+1 vid(end_flag0+1)'-vid(start_flag0)' ...
        zeros(length(start_flag0),1) i*ones(length(start_flag0),1)]];
end
if ~isempty(start_flag1)
    Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [ones(length(start_flag1),1) ...
        start_flag1' end_flag1' end_flag1'-start_flag1'+1 vid(end_flag1+1)'-vid(start_flag1)' ...
        zeros(length(start_flag1),1) i*ones(length(start_flag1),1)]];
end
if ~isempty(start_flag2) & start_flag2~=start_flag_2
    Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [2 start_flag2 end_flag2 ...
        end_flag2-start_flag2+1 vid(end_flag2+1)-vid(start_flag2) 0 i]];
end
if ~isempty(start_flag3) & start_flag3~=start_flag_3
    Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [3 start_flag3 end_flag3 ...
        end_flag3-start_flag3+1 vid(end_flag3+1)-vid(start_flag3) 0 i]];
end
end
    
Intervals.CumulativeIntervalArray=sortrows(Intervals.CumulativeIntervalArray,[7,2]);

cia=Intervals.CumulativeIntervalArray;
logik=[];
logik=(cia(:,1)==-2);
Time2FirstEvent=cia(logik,5); 
time_empty=max(Time2FirstEvent)*ones(length(Time2FirstEvent),1);
t=(Time2FirstEvent~=time_empty);
time21st=Time2FirstEvent(t);
mean(time21st)


cia=Intervals.CumulativeIntervalArray;
logik=[];
logik=cia(:,1)==1; % &cia(:,4)>1;
DwellTime=cia(logik,5);
mean(DwellTime)

% pd = fitdist(DwellTime,'Exponential')
% x_values = 0:1:250;
% y = pdf(pd,x_values);
% plot(x_values,y)
% hold on
% histogram(DwellTime,round(size(DwellTime,1)/100),'Normalization','probability')

%save IDS2.mat Intervals;

logik1st=cia(:,1)==-2; % this is to sort in order of first arrival time
Arrivals=[cia(logik1st,3) cia(logik1st,7)]; %1st column is arrival frame, 2nd column is aoi#
Sorted_1stevent=sortrows(Arrivals,1); % Sort in the order of  arrival time for 1st arrival
newcia=[]; % initialize new cumulative interval array for sorted IDS
for i=1:max(size(Sorted_1stevent));
    logik=cia(:,7)==Sorted_1stevent(i,2);
    instant_cia=cia(logik,:); instant_cia(:,7)=i;
    newcia=[newcia;instant_cia]; 
end
Intervals.CumulativeIntervalArray=newcia;

% imscroll_off=load('DwellTime_imscroll_test_2016_0127.mat', '-mat');
% imscroll_on=load('time21st_imscroll_test_2016_0127.mat', '-mat');
 
cdfplot(DwellTime)
hold on
% cdfplot(imscroll_off.DwellTime)
xlabel('seconds');
legend('pipeline DwellTime', 'imscroll DwellTime', 'Location', 'SouthEast')
hold off
figure;
cdfplot(time21st)
hold on
% cdfplot(imscroll_on.time21st)
xlabel('seconds');
legend('pipeline time to 1st event', 'imscroll time to 1st event', 'Location', 'SouthEast')
hold off

num_traces=max(Intervals.CumulativeIntervalArray(:,7));
save IDS2_sorted.mat Intervals;
vid.ttb = (1:size(traces2_allfr,2))*1000;
save vidfile.mat vid
% (file_path,header_path,min_frame,max_frame,spot_number,nodata,nopeak,peak,silent)
func_graph_intervals_nooverlap_v1pt0('IDS2_sorted','vidfile.mat',1,size(traces2_allfr,2),num_traces,[1 1 1],[1 1 1],[0 0 0],0);
% func_graph_intervals_nooverlap_v1pt0('sorted_1301_IDS_exp003_561nm_1rad_3x6Up_5pxls_140low_250high','vidfile.mat',1,599,500,[1 1 1],[1 1 1],[0 0 0],0);

% [fn fp]=uigetfile
% eval(['load ' [fp fn] ' -mat'])
cia=Intervals.CumulativeIntervalArray;
logik=[];
logik=cia(:,1)==1;
DwellTime=cia(logik,5);
% The following removes 1-frame intervals
logik1=[];
logik1=DwellTime(:)~=1;
DwellTimeNo1fr=DwellTime(logik1);
% figure(32); optN = sshist(DwellTimeNo1fr); hist(DwellTimeNo1fr,optN);

% to make a cdf from dwell times
start_time=0; % these is start and  
end_time=5; % end time in seconds.
upper=max(DwellTime);
cdf=[];
for i=start_time:0.1:end_time
  A=DwellTime<=i;
  B=sum(A);
  cdf=[cdf;[i B]];

end

t=cdf(:,1);
risc=cdf(:,2);
dark=cdf(:,2);
dfittool(DwellTime)

%%

% 
% 
% %%
% 
% plot(cumsum(sv'),'linewidth',2)
% legend('x-drift','y-drift')
% 
% ylabel('Cumultive Drift [pixel]')
% xlabel('t [frames]')
% % grid on
% 
%  y = cumsum(sv');
%  
% export_fig drifNG.eps -transparent 
% 
% 
% export_fig drifNG.png -transparent 
% %%
% a = histogram(rawFitResultsCam1.Photons,15);
% 
% ylabel('Frequency [#]')
% xlabel('Photons')
% grid on
% export_fig hist.eps -transparent 
% export_fig hist.png -transparent 



%%

DetectedCoords = rawInitialFitResultsCam1.Coord; %(rawInitialFitResultsCam1.Frame ==0,:);
tDetectedCoords = transformPointsInverse(tformTotal,DetectedCoords);
idx = C(1,1) < tDetectedCoords(:,1) &...
C(1,1)+C(2,1) > tDetectedCoords(:,1)&...
C(1,2) < tDetectedCoords(:,2) &...
C(1,2)+C(2,2) > tDetectedCoords(:,2);

dipshow(dataCam1(:,:,1))
% hold on
% plot(tDetectedCoords(idx,1),tDetectedCoords(idx,2),'xg')
hold on
plot(DetectedCoords(idx,1),DetectedCoords(idx,2),'xr')

DetectedCoords = DetectedCoords(idx,:);

gridcell=10; % This is the size of the dark AOI grid cell (smaller values lead to 

% higher grid density of dark AOIs; this should be greater or equal to the AOI size
mindist=10; % This is the minimal distance in pixels from non-dark AOIs to maintain 

xmin=min(DetectedCoords(:,1))+gridcell/2;
xmax=max(DetectedCoords(:,1))-gridcell/2;
ymin=min(DetectedCoords(:,2))+gridcell/2;
ymax=max(DetectedCoords(:,2))-gridcell/2;
xgrid=linspace(xmin,xmax,(xmax-xmin)/gridcell);
ygrid=linspace(ymin,ymax,(ymax-ymin)/gridcell);

[x,y]=meshgrid(xgrid,ygrid);
a=ones(size(x));

for i=1:length(DetectedCoords(:,1)) 
    xnotdark=find(abs(DetectedCoords(i,1)-xgrid)<mindist);
    ynotdark=find(abs(DetectedCoords(i,2)-ygrid)<mindist);
    a(ynotdark,xnotdark)=0;
end
darkCoords=[];
% Draw dark locations (blue) and RNA target locations (red)
darkCoords(:,1)=x(a==1);
darkCoords(:,2)=y(a==1);
disp([num2str(numel(darkCoords(:,1))) ' dark AOIs found'])

dipshow(dataCam1(:,:,1))
hold on
scatter(darkCoords(:,1),darkCoords(:,2), 'xb')
scatter(DetectedCoords(:,1),DetectedCoords(:,2),'xr')

dipshow(dataCam2(:,:,1))
hold on
tdarkCoords = transformPointsInverse(tformTotal,darkCoords);
scatter(tdarkCoords(:,1),tdarkCoords(:,2), 'xb')
scatter(tDetectedCoords(idx,1),tDetectedCoords(idx,2),'xr')

coordsDarkCam1=[];
coordsDarkCam2=[];
idxLog=[];
for w=1:size(dataCam1,3)
    coords_dark_nth = darkCoords-repmat(drift(w,:),size(darkCoords,1),1);
    coordsDarkCam1=cat(1,coordsDarkCam1,[coords_dark_nth (w-1).*ones(size(darkCoords,1),1)]);
    coordsDarkCam2=cat(1,coordsDarkCam2,[transformPointsInverse(tformTotal,coords_dark_nth) (w-1).*ones(size(darkCoords,1),1)]); 
    idxLogCam1(:,:,w) = logical(maxf(coord2image(coords_dark_nth,[size(dataCam1,1),size(dataCam1,2)]),6,'rectangular'));
    idxLogCam2(:,:,w) = logical(maxf(coord2image(transformPointsInverse(tformTotal,coords_dark_nth) ,[size(dataCam1,1),size(dataCam1,2)]),6,'rectangular'));
end

% [canCoordsCam1,detParCam1,cutProcessCam1] = LLRMapv2(dataCam1,PSFSigma,[],idxLogCam1,0.05,8,true);
[canCoordsCam2,detParCam2,cutProcessCam1] = LLRMapv2(dataCam2,PSFSigma,[],idxLogCam2,0.05,8,true);


%%
canCoordsCam2(:,3) = ones(size(canCoordsCam2,1),1);
img = coord2image(canCoordsCam2(:,1:2),[size(cutProcessCam1,1),size(cutProcessCam1,2)]);
msr = measure(squeeze(label(logical(img))), ones(size(img)), ({'Gravity'}));
coords =msr.Gravity';

[image_out,R] = affine_trans(dataCam1(:,:,1), out(1:2), out(3:4), out(5));
dipshow(image_out)
hold on
canCoordsCam2UnCut=round(canCoordsCam2+(1.5*(2*PSFSigma+1)-0.5).*[ones(size(canCoordsCam2,1),1) ones(size(canCoordsCam2,1),1) zeros(size(canCoordsCam2,1),1)]);
scatter(tdarkCoords(:,1),tdarkCoords(:,2), 'o')
plot(canCoordsCam2UnCut(:,1),canCoordsCam2UnCut(:,2),'x')
canCoordsCam2UnCutSum=round(coords+(1.5*(2*PSFSigma+1)-0.5));

plot(canCoordsCam2UnCutSum(:,1),canCoordsCam2UnCutSum(:,2),'d')

paramsPreFilterFits = getDefaultParamsPreFilterFits;
paramsPreFilterFits.minPixelDist=8;
[ maskPreFiltCam1 ] =  preFilterFits([ canCoordsCam2UnCutSum ones(size(canCoordsCam2UnCutSum,1),1)] ,[],paramsPreFilterFits);
canCoordsCam1UnCutSum = transformPointsForward(tformTotal,canCoordsCam2UnCutSum);
canCoordsCam1UnCutSum=canCoordsCam1UnCutSum(maskPreFiltCam1,:);
canCoordsCam2UnCutSum=canCoordsCam2UnCutSum(maskPreFiltCam1,:);

%%
dipshow(dataCam1(:,:,1))
hold on
plot(canCoordsCam1UnCutSum(:,1),canCoordsCam1UnCutSum(:,2),'d')

paramsFit = getDefaultParamsFit;
paramsFit.Iterations=16;
paramsFit.FitSigma=true;

canCoordsCam1Dark=[];
canCoordsCam2Dark=[];
for w=1:size(dataCam1,3)
    coords_nth = canCoordsCam1UnCutSum-repmat(drift(w,:),size(canCoordsCam1UnCutSum,1),1);
    canCoordsCam1Dark=cat(1,canCoordsCam1Dark,[[coords_nth(:,1) coords_nth(:,2)] (w-1).*ones(size(coords_nth,1),1)]);
    canCoordsCam2Dark=cat(1,canCoordsCam2Dark,[transformPointsInverse(tformTotal,[coords_nth(:,1) coords_nth(:,2)]) (w-1).*ones(size(coords_nth,1),1)]); 
end

x = reshape(canCoordsCam2Dark(:,1),[33 600]);
y = reshape(canCoordsCam2Dark(:,2),[33 600]);
t = reshape(canCoordsCam2Dark(:,3),[33 600]);


x1 = reshape(canCoordsCam1Dark(:,1),[33 600]);
y1 = reshape(canCoordsCam1Dark(:,2),[33 600]);
t1 = reshape(canCoordsCam1Dark(:,3),[33 600]);
clear rawFitResultsDarkCam1
clear rawFitResultsDarkCam2

shp = randperm(33,33);
[ rawFitResultsDarkCam1 ] = fitBoxCenters( single(squeeze(dataCam1)),[reshape(x1(shp,:),[numel(shp)*size(x1,2) 1]) reshape(y1(shp,:),[numel(shp)*size(y1,2) 1]) reshape(t1(shp,:),[numel(shp)*size(t1,2) 1])],paramsFit); %canCoordsCamshpDark
[ rawFitResultsDarkCam2 ] = fitBoxCenters( single(squeeze(dataCam2)),[reshape(x(shp,:),[numel(shp)*size(x,2) 1]) reshape(y(shp,:),[numel(shp)*size(y,2) 1]) reshape(t(shp,:),[numel(shp)*size(t,2) 1])],paramsFit);

dipshow(rawFitResultsDarkCam2.ROIStack)

optionsLLR = dipSubLoc2DSetOptions;
optionsLLR.BoxCenters = [rawFitResultsDarkCam1.Coord(:,1) rawFitResultsDarkCam1.Coord(:,2) rawFitResultsDarkCam1.Frame(:)];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 
index = double(floor(255.*(detParCam1.PH1(:,1)))+1);
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = dataCam1; %dip_image(stretch(permute( single(squeeze(stretch(cut(dataCam1,[size(cutProcessCam1,2) size(cutProcessCam1,2) size(cutProcessCam1,3)])))),[1 2 3])),'uint8')';
h = dipSubLoc2D(optionsLLR);  
hold on
plot(canCoordsCam1UnCutSum(:,1),canCoordsCam1UnCutSum(:,2),'d')



optionsLLR = dipSubLoc2DSetOptions;
optionsLLR.BoxCenters = [rawFitResultsDarkCam2.Coord(:,1) rawFitResultsDarkCam2.Coord(:,2) rawFitResultsDarkCam2.Frame(:)];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 
index = double(floor(255.*(ones(size(detParCam1.PH1(:,1))))+1));
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = dataCam2; %dip_image(stretch(permute( single(squeeze(stretch(cut(dataCam1,[size(cutProcessCam1,2) size(cutProcessCam1,2) size(cutProcessCam1,3)])))),[1 2 3])),'uint8')';
h = dipSubLoc2D(optionsLLR);  
hold on
plot(canCoordsCam2UnCutSum(:,1),canCoordsCam2UnCutSum(:,2),'d')



[traces2_allfr_dark,seqDark,  transmat1, prior1, mu1, Sigma1, mixmat1] =  HMMKEstimate(rawFitResultsDarkCam1,rawFitResultsDarkCam2,[reshape(x(shp,:),[numel(shp)*size(x,2) 1]) reshape(y(shp,:),[numel(shp)*size(y,2) 1]) reshape(t(shp,:),[numel(shp)*size(t,2) 1])],tformTotal);

makeMovieTraces(traces2_allfr_dark,seqDark, rawFitResultsDarkCam1,rawFitResultsDarkCam2,canCoordsCam2Dark, movieFileName2,100)
%%



NFrames  =rawFitResultsDarkCam2.Frame(end)+1;
maskFilt2 = true(size(rawFitResultsDarkCam2.Photons,1),1);
maskFilt1 = true(size(rawFitResultsDarkCam1.Photons,1),1);

photons1 = rawFitResultsDarkCam1.Photons(maskFilt1,:);
% frames2 = rawFitResultsDarkCam2.Frame(maskFilt2,:);
% frames1 = rawFitResultsDarkCam1.Frame(maskFilt1,:);
bg2 = rawFitResultsDarkCam2.Bg(maskFilt2,:);
photons2 = rawFitResultsDarkCam2.Photons(maskFilt2,:);
coords1 = rawFitResultsDarkCam2.Coord(maskFilt2,:);
coords2 = rawFitResultsDarkCam2.Coord(maskFilt2,:);

photons2I = reshape(photons2,[size(photons1,1)/NFrames NFrames]);
bg2I = reshape(bg2,[size(photons2,1)/NFrames NFrames]);

photons1I = reshape(photons1,[size(photons1,1)/NFrames NFrames]);
coords1x = reshape(coords1(:,1),[size(photons1,1)/NFrames NFrames]);
coords1y = reshape(coords1(:,2),[size(photons1,1)/NFrames NFrames]);


coords2x = reshape(coords2(:,1),[size(photons1,1)/NFrames NFrames]);
coords2y = reshape(coords2(:,2),[size(photons1,1)/NFrames NFrames]);

esigmax = reshape(rawFitResultsDarkCam2.Sigma(:,1),[size(photons1,1)/NFrames NFrames]);
esigma = reshape(rawFitResultsDarkCam2.Sigma(:,2),[size(photons1,1)/NFrames NFrames]);

% coordsCam2x = reshape(coordsCam2(:,1),[size(photons1,1)/NFrames NFrames]);
% coordsCam2y = reshape(coordsCam2(:,2),[size(photons1,1)/NFrames NFrames]);


%%

% delta = transformPointsInverse(tformTotal,rawFitResultsDarkCam1.Coord)-rawFitResultsDarkCam2.Coord;
% deltax = reshape(delta(:,1),[size(delta(:,1),1)/NFrames NFrames]);
% deltay = reshape(delta(:,2),[size(delta(:,2),1)/NFrames NFrames]);
PSFSigma=1.39;
seq=[];
% seq(:,:,1) = sqrt((coords2x-coordsCam2x).^2+ (coords2y-coordsCam2y ).^2) ;
% seq(:,:,5) = min(ecrlb,10);
seq(:,:,1) = min(esigma,10);
seq(:,:,2) = min(esigmax,10);
seq(:,:,3) = photons2I-bg2I*(PSFSigma*2)^2;
% seq(:,:,4) = coords2x-coordsCam2x ; %deltax;
% seq(:,:,4)=seq(:,:,4)-repmat(mean(seq(:,:,4),2),[1 size(seq,2)]);
% seq(:,:,5) = coords2y-coordsCam2y ; %deltay; 
% seq(:,:,5) =seq(:,:,5)-repmat(mean(seq(:,:,5),2),[1 size(seq,2)]);

O = size(seq,3); %opservation states
T = size(seq,2); %time
nex = size(seq,1); %tracks
M = 1; %g muxture of output
Q = 2; %states

data = permute(seq,[3 2 1]);
% data = randn(O,T,nex);

% initial guess of parameters
prior0 = normalise(rand(Q,1));
transmat0 = mk_stochastic(rand(Q,Q));

cov_type = 'full';
idx = (photons2I > mean(mean(bg2I*(PSFSigma*2)^2)) +std(std(bg2I*(PSFSigma*2)^2)));

[mu0, Sigma0] = mixgauss_init(Q*M, data, cov_type);
% mu0
% 
% 
% mu0(1,1) = mean(min(esigma(idx),10));  
% mu0(2,1) = mean(min(esigmax(idx),10)); 
% mu0(3,1) = mean(photons2I(idx)); 
% mu0(1,2) = mean(min(esigma(~idx),10));  
% mu0(2,2) = mean(min(esigmax(~idx),10)); 
% mu0(3,2) = mean(photons2I(~idx)); 
% 
% mu0
% 
% % 
% % mu0
% [v,I] = sort(mu0(3,:))
% mu0=mu0(:,I);
% Sigma0=Sigma0(:,:,I);
% mu0 = reshape(mu0, [O Q M]);
% 
% 
% Sigma(:,:,2) = cov([min(esigma(~idx),10) min(esigmax(~idx),10) photons2I(~idx)]);
% Sigma(:,:,1) = cov([min(esigma(idx),10) min(esigmax(idx),10) photons2I(idx)]);
% 
% if isposdef(Sigma(:,:,1))
%     Sigma0(:,:,1) = Sigma(:,:,1);
% end
% if isposdef(Sigma(:,:,2))
%     Sigma0(:,:,2) = Sigma(:,:,2);
% end

% 

Sigma0 = reshape(Sigma0, [O O Q M]);
mixmat0 = mk_stochastic(rand(Q,M));

% 'adj_prior',0
% 'adj_trans' - if 0, do not change transmat [1]
% 'adj_mix' - if 0, do not change mixmat [1]
% 'adj_mu' - if 0, do not change mu [1]
% 'adj_Sigma' - if 0, do not change Sigma [1]
onlyTrans = true ;
if onlyTrans
    [LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
    mhmm_em(data, prior1, transmat1,  mu0, Sigma0, mixmat1, 'max_iter', 10,'adj_mu',0,'adj_Sigma',0);
else
[LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
    mhmm_em(data, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 10);
end

loglik = mhmm_logprob(data, prior1, transmat1, mu1, Sigma1, mixmat1);

for ex=1:size(data,3)
    B = mixgauss_prob(data(:,:,ex), mu1, Sigma1, mixmat1);
    [path(ex,:)] = viterbi_path(prior1, transmat1, B);
end
traces2_allfr_dark=(path-1);


makeMovieTraces(traces2_allfr_dark,seq, rawFitResultsDarkCam1,rawFitResultsDarkCam2,canCoordsCam2Dark, movieFileName2,100)


