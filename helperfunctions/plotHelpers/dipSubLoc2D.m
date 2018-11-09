function varargout = dipSubLoc2D(varargin)
% dipSubLoc2D     plot subregions and localizations
%
% USAGE:
%   h = dipSubLoc2D(options); %setup dipSubLoc2D figure
%   h = dipSubLoc2D(fh); %setup listener for dipSubLoc2D figure



if nargin > 2
    error('dipSubLoc2D:ToManyInputs','dipSubLoc2D: 0 to 2 inputs required');
end

switch nargin
    case 0
        options = dipSubLoc2DSetOptions;
    case 1
        if ishandle(varargin{1})
            h = varargin{1};
            udata = get(h,'userdata');
            if isfield(udata,'dipSubLoc2DData')
                udata.dipSubLoc2DData.h = h;
                addlistener(udata.dipSubLoc2DData.h,'UserData','PostSet',@curslice);
                set(h,'userdata',udata);
                return;
            else
                error('dipSubLoc2D:InputMustBeDipTrack','dipSubLoc2D: input figure must be initialized using dipSubLoc2D')
            end
        end
        options = varargin{1};
end

%initilize figure
if isempty(options.h)
    h = figure;
else
    h = options.h;
end
%show image
if isempty(options.im)
    warning('dipSubLoc2D:NoIm','dipSubLoc2D: options.im is empty initializing with newim(256,256,10)')
    if ~isempty(options.BoxCenters)
        dipshow(h,newim(256,256,max(options.BoxCenters(:,3))));
    else
        dipshow(h,newim(256,256,10));
    end
else
    if ishandle(options.h) && strcmp(get(options.h,'type'),'figure')       
        dipshow(h,options.im);
    elseif ishandle(options.h)
        imshow(mat2gray(options.im))
    end
end
%initialize dipData
dipSubLoc2DData.BoxCenters = options.BoxCenters;
dipSubLoc2DData.BoxSize = options.BoxSize;
dipSubLoc2DData.BoxColor = options.BoxColor;
dipSubLoc2DData.plotBoxes = options.plotBoxes;
dipSubLoc2DData.fitCoords = options.fitCoords;
dipSubLoc2DData.fitCoordsMarker = options.fitCoordsMarker;
dipSubLoc2DData.fitColor = options.fitColor;
dipSubLoc2DData.markersize = options.markersize;
dipSubLoc2DData.linewidth = options.linewidth;
dipSubLoc2DData.h = h;
dipSubLoc2DData.ha = findall(h,'type','axes');
dipSubLoc2DData.dipTrackObj = dipTrackObj(0);
%add dipSubLoc2DData to userdata
udata = get(h,'userdata');
udata.dipSubLoc2DData = dipSubLoc2DData;
set(h,'userdata',udata);
%add listener
addlistener(udata.dipSubLoc2DData.h,'UserData','PostSet',@curslice);

dipmapping(h,'slice',1);
dipmapping(h,'slice',0);

%output figure handle
if nargout == 1
    varargout{1} = h;
end


function update_boxes(udata)
%update boxes for subregions in dipSubLoc2D figure

%don't update if slice is the same
if isempty(udata.slicing) || isempty(udata.dipSubLoc2DData.BoxCenters)
    return;
end
%find boxes in current slice
boxIdx = find(udata.dipSubLoc2DData.BoxCenters(:,3) == udata.curslice);
%random color scheme
if isempty(udata.dipSubLoc2DData.BoxColor)
    cTmp = jet(length(boxIdx));
    [v idx] = sort(rand(size(boxIdx)));
    c = cTmp(idx,:);
else if size(udata.dipSubLoc2DData.BoxColor,1)
        c = repmat(udata.dipSubLoc2DData.BoxColor,[length(boxIdx) 1]);
    else
        c = udata.dipSubLoc2DData.BoxColor(boxIdx,:);
    end
end
%delete all current boxes
delete(findall(udata.dipSubLoc2DData.ha,'tag','plotBox'));
delete(findall(udata.dipSubLoc2DData.ha,'tag','dipTrackLocMarker'));
if udata.dipSubLoc2DData.plotBoxes && ~isempty(udata.dipSubLoc2DData.BoxCenters)
    %plot boxes
    plotBox(udata.dipSubLoc2DData.BoxCenters(boxIdx,1),udata.dipSubLoc2DData.BoxCenters(boxIdx,2),...
        c,udata.dipSubLoc2DData.BoxSize(1),udata.dipSubLoc2DData.BoxSize(2),udata.dipSubLoc2DData.linewidth,...
        udata.dipSubLoc2DData.ha);
end
if udata.dipSubLoc2DData.fitCoordsMarker && ~isempty(udata.dipSubLoc2DData.fitCoords)
    if isempty(udata.dipSubLoc2DData.fitColor)
        for ii = 1:length(boxIdx)
            fitIdx = find(udata.dipSubLoc2DData.fitCoords(:,3) == boxIdx(ii));
            for jj = 1:length(fitIdx)
                %plot fit coordinates
                plotLocMarker(udata.dipSubLoc2DData.fitCoords(fitIdx(jj),1),udata.dipSubLoc2DData.fitCoords(fitIdx(jj),2),...
                    c(ii,:),udata.dipSubLoc2DData.markersize,udata.dipSubLoc2DData.ha)
            end
        end
    else
       for ii = 1:length(boxIdx)
            fitIdx = find(udata.dipSubLoc2DData.fitCoords(:,3) == boxIdx(ii));
            c = udata.dipSubLoc2DData.fitColor(fitIdx,:);
            for jj = 1:length(fitIdx)
                %plot fit coordinates
                plotLocMarker(udata.dipSubLoc2DData.fitCoords(fitIdx(jj),1),udata.dipSubLoc2DData.fitCoords(fitIdx(jj),2),...
                    c(jj,:),udata.dipSubLoc2DData.markersize,udata.dipSubLoc2DData.ha)
            end
       end
    end
end


function curslice(h,evnt)
%check updating of curslice field in userdata 

%get userdata

if isobject(evnt)
    udata = get(evnt.AffectedObject,'UserData');
else
    udata = get(evnt,'NewValue');
end

%get current slice
if isfield(udata,'curslice') && udata.dipSubLoc2DData.dipTrackObj.slice ~= udata.curslice
    udata.dipSubLoc2DData.dipTrackObj.slice = udata.curslice; %update slice property
    update_boxes(udata); %update figure
    %     set(udata.dipTrackData.h,'userdata',udata)
end


function h = plotBox(x,y,c,w,h,lw,ha)
% PLOTBOX   add box around a given coordinate
% 
% h = plotBox(x,y,c,w,h,lw)
% 
% INPUTS
%   x - x coordinate(s)
%   y - y coordinate(s)
%   c - color of box
%   w - width of box
%   h - height of box
%   lw - linewidth
%   ha - axes handle
% OUTPUT
%   h - handle for line object
% 
% Created by Pat Cutler October 2009


if size(c,1) == 1
    c = repmat(c,[size(x,1) 1]);
end
distx = w/2;
disty = h/2;
if size(x,1) == 1
    x = shiftdim(x);
    y = shiftdim(y);
end
xc = [x-distx x+distx x+distx x-distx x-distx];
yc = [y+disty y+disty y-disty y-disty y+disty];
% h = line(xc,yc,'color',c,'linewidth',lw,'parent',ha);

hold on
for ii = 1:length(x)
    h(ii) = line(xc(ii,:),yc(ii,:),'color',c(ii,:),'linewidth',lw,'parent',ha,'tag','plotBox');
end
hold off


function plotLocMarker(x,y,c,markersize,ha)
%plot localization marker

line(x,y,...
    'marker','o','color',c*0.8,'markersize',markersize(1),...
    'tag','dipTrackLocMarker','parent',ha);
line(x,y,...
    'marker','o','color',c*0.6,'markersize',markersize(2),...
    'markerfacecolor',c*.3,'tag','dipTrackLocMarker','parent',ha);

% hold(ha,'on')
% scatter(ha,x,y,...
%     markersize(1),c,'filled','marker','o',...
%     'markeredgecolor',[1 1 1],'tag','dipTrackLocMarker');
% scatter(ha,x,y,...
%     markersize(2),c,'filled','marker','o',...
%     'markeredgecolor',[0 0 0],'tag','dipTrackLocMarker');
% hold(ha,'off')