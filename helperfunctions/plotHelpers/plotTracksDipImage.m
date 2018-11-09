function [fh fhcb] = plotTracksDipImage(obj,colorBarFig,saveStr,trackNum,tracksVal,minMaxVal,cmap,valueName)
% PLOTTRACKSDIPIMAGE  make dipTrack figure
%
% USAGE:
%   obj.plotTracksDipImage;
%   obj.plotTracksDipImage(obj,colorBarFig,saveStr,tracksNum,tracksNum,tracksVal,minMaxVal);
%
% INPUTS:
%   colorBarFig: (binary) 1 make colorbar figure.
%   saveStr: string with appendix for saveName. if empty don't save. Default empty.
%   trackNum: (array) with same dimensions as obj.Connect.tracks
%           if not specified then obj.Connect.tracks used.
%   tracksVal: (array) with same dimensions as first 2 dimentsion of obj.Connect.tracks
%             if not specified then track number is used.
%   minMaxVal: 1 by 2 vector with min max values for color
%              scaling
%   cmap: color map
%   valueName: (string) describing the values in trackVal.
%                Default is 'track#'.


%get tracks
[tracks, ~, subRegionModel] = obj.getTracks;
% addpath('C:\Users\smithc5\Dropbox\grunwaldLabPc\summary\detection\backup\helperfunctions')
% fileName{1} = 'LoGforCarlasTH5000.trk';
% pathName{1} = 'C:\Users\smithc5\Dropbox\grunwaldLabPc\summary\detection\backup\localize\';
% [ ~,TracksRaw ] = importTrackFile(fileName,pathName, 1 )
% TracksRaw{1}.t

% tracks = zeros(numel(obj.Tracks),obj.Stats.Data.size(3),ndims);
% TracksRaw
% subRegionModel = zeros(numel(obj.Tracks),obj.Stats.Data.size(3),4);
%check inputs
if ~exist('colorBarFig','var') || isempty(colorBarFig)
    colorBarFig = 0;
end
if ~exist('saveStr','var')
    saveStr = [];
end
if ~exist('trackNum','var') || isempty(trackNum)
    trackNum = 1:size(tracks,1);
end
if ~exist('tracksVal','var') || isempty(tracksVal)
    tracksVal = logical(tracks(:,:,1)).*repmat((1:size(tracks,1))',[1 size(tracks,2)]);
end
if ~exist('valueName','var') || isempty(valueName)
    valueName = 'track#';
end
if ~exist('minMaxVal','var')
    minMaxVal = [];
end

%get/make figure
if isfield(obj.FigInfo,'plotTracksDipImage') && ~isempty(obj.FigInfo.plotTracksDipImage)...
        && isfield(obj.FigInfo.plotTracksDipImage,'fh') && ~isempty(obj.FigInfo.plotTracksDipImage.fh)...
        && ishandle(obj.FigInfo.plotTracksDipImage.fh)
    fh = obj.FigInfo.plotTracksDipImage.fh;
    clf(fh)
    set(fh,'userdata',[],'tag','')
else
    fh = figure;
end
if colorBarFig
    if isfield(obj.FigInfo,'plotTracksDipImage') && ~isempty(obj.FigInfo.plotTracksDipImage)...
            && isfield(obj.FigInfo.plotTracksDipImage,'fhcb') && ~isempty(obj.FigInfo.plotTracksDipImage.fhcb)...
            && ishandle(obj.FigInfo.plotTracksDipImage.fhcb)
        fhcb = obj.FigInfo.plotTracksDipImage.fhcb;
    else
        fhcb = figure;
    end
else
    fhcb = [];
end
%setup dipTrack
dipTrackOptions = dipTrackSetOptions;
dipTrackOptions.h = fh;
if ~isempty(strfind(get(fh,'tag'),'DIP_Image'))
    dipTrackOptions.updateImage = 0;
else
    if isempty(obj.Data)
        obj.loadData;
    end
    dipTrackOptions.im = uint8(stretch(obj.Data,0,100,0,255));
end

udata = get(fh,'userdata');
if isfield(udata,'dipTrackData');
    dipTrackOptions.plotTracksOptions = udata.dipTrackData.plotTracksOptions;
    nplotTracks = numel(dipTrackOptions.plotTracksOptions)+1;
    dipTrackOptions.plotTracksOptions(nplotTracks) = plotTracksV1SetOptions;
else
    nplotTracks = 1;
end

dipTrackOptions.plotTracksOptions(nplotTracks).trackNum = trackNum;
dipTrackOptions.plotTracksOptions(nplotTracks).tracks = tracks;
dipTrackOptions.plotTracksOptions(nplotTracks).subRegionModel = subRegionModel;
dipTrackOptions.plotTracksOptions(nplotTracks).tracksVal = tracksVal;
dipTrackOptions.plotTracksOptions(nplotTracks).minMaxVal = minMaxVal;
dipTrackOptions.plotTracksOptions(nplotTracks).marker = 'none';
dipTrackOptions.plotTracksOptions(nplotTracks).addShadow = dipTrackOptions.plotTracksOptions.linewidth;

if exist('cmap','var') && ~isempty(cmap)
    dipTrackOptions.plotTracksOptions(nplotTracks).cmap = cmap;
else
    dipTrackOptions.plotTracksOptions(nplotTracks).cmap = lines(size(tracks,1));
end
dipTrackOptions.h = fh;
dipTrackOptions.headMarkersize = [6 3];
% dipTrackOptions.npoints = inf;
%make figure
obj.FigInfo.plotTracksDipImage.fh = dipTrack(dipTrackOptions);
% %change size
% diptruesize(fh,200);
%ensure that figure will fit on monitor at 1000%
monPos = get(0,'MonitorPositions');
tmpSz = obj.Stats.Data.size([1 2])*10;
if any(tmpSz([1 2]) > min(monPos(:,[3 4])))
    diptruesize(fh,min(floor(min(monPos(:,[3 4]))./min(tmpSz([1 2]))*10-1))*100);
else
    %change size
    diptruesize(fh,1000);
end
%add scaleLine
delete(findall(fh,'tag','scaleLine'))
scaleLine(findall(fh,'type','axes'),2,obj.Stats.Data.PixelSize);

if colorBarFig
    figure(fhcb)
    hcb = plot_colorbar([100 20],'v','',jet,[1 size(tracks,1)]);
    title(valueName)
    set(findall(gca,'type','text'),'fontsize',16,'fontweight','bold')
    set(gca,'fontsize',16,'fontweight','bold')
end

if ~isempty(saveStr) %save figure
    obj.saveFig(fh,'plotTracksDipImage',saveStr);
    if colorBarFig
        obj.saveFig(fhcb,'plotTracksDipImage',[saveStr '-colorbar']);
    end
end
