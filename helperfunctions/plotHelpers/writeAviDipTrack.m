function writeAviDipTrack(movieFileName, tracks,sequence, scaleLineLength,pixelSize,TimeStep, xzoom, yzoom,tzoom)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Optional
% xzoom = []; %make empty for full image
% yzoom = []; %make empty for full image
% vidObj.FrameRate = round(1/TimeStep);

%% example for making an avi
% pixelSize = 1;
% TimeStep=1;
% path = 'C:\Users\smithc5\Dropbox\grunwaldLabPc\summary\detection\'; 
% movieFileName = fullfile(path,'dipTrackTstMovie.avi');
% addpath('C:\Users\smithc5\Dropbox\grunwaldLabPc\summary\detection\backup\helperfunctions')
% fileName{1} = 'LoGforCarlasTH5000.trk';
% pathName{1} = 'C:\Users\smithc5\Dropbox\grunwaldLabPc\summary\detection\backup\localize\';
% [ ~,~,tracks ] = importTrackFile(fileName,pathName, 1 )

%length of scaleLine in microns
% scaleLineLength = 1; 

if nargin < 7 || isempty(xzoom) || isempty(yzoom) || isempty(tzoom)
    xzoom = []; %make empty for full image
    yzoom = []; %make empty for full image
    tzoom=size(tracks,2);
end

%set limits for ROI
optoins = plotTracksV1SetOptions('colorbar',0,'tracks',tracks,'tracksVal',repmat([ceil(size(tracks,2).*rand(1,size(tracks,1)))]',[1 , size(tracks,2)]));
 plotTracksV1(optoins)

dipTrackOptions = dipTrackSetOptions;
dipTrackOptions.plotTracksOptions = optoins;
dipTrackOptions.im = sequence;
dipTrackOptions.plotTracksOptions.colorLineByTrackNum = 1;
dipTrack(dipTrackOptions)
h=dipTrack(dipTrackOptions);

%resize figure
diptruesize(h,300)

%zoom into ROI
pause(0.1)
ha = findall(h,'type','axes');
%set xlimits
if ~isempty(xzoom)
    xlim(ha,xzoom+[-0.5 0.5])
end
if ~isempty(yzoom)
    ylim(ha,yzoom+[-0.5 0.5])
end

%add appropriate scaleLine
delete(findall(h,'tag','scaleLine')) %delete current scaleLine
% scaleLine(ha,scaleLineLength,pixelSize);

clear mex
% Prepare VideoWriter object
vidObj = VideoWriter(movieFileName);
vidObj.FrameRate = round(TimeStep);
%open avi file
open(vidObj);

for ii = 1:min(size(tracks,2),tzoom)
    %update frame
    dipmapping(h,'slice',ii-1)
    % Write each frame to the file.
    currFrame = getframe(ha);
    data = uint8(extend(currFrame.cdata,[194 194 3]));
    currFrame.cdata = data;
    writeVideo(vidObj,currFrame);
end

%close avi file
close(vidObj);
close(h);
end

