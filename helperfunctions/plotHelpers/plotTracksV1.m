function [h ha hp] = plotTracksV1(options)
% PLOTTRACKSV1    plot tracks using specified options
%
% h = plotTracksV1(options)
%
% INPUTS
%   options - options structure with input options. See
%             plotTracksSetOptions for more details.
% OUTPUTS
%   h - figure handle
%   ha - axes handle
%   hp - patch handle
%


%%
% options = plotTracksSetOptions;
% options.tracks = tracks;
% options.tracksVal = tracksWv;
% options.t = 1;
%%

% highlight tracks if input is str
if ischar(options)
    try
        eval(options);
        return
    catch err
        getReport(err)
        error('plotTracks:ImproperInput','plotTracks: character input is not correct')
    end
    %     if strfind(options,'highlightTracks')
    %         h = eval(options);
    %         return;
    %     else if strfind(options,'initDataTip')
    %             eval(options);
    %             return;
    %         else if strfind(options,'initDataTip')
    %                 eval(options,'resetTracks')
    %             else
    %                 error('plotTracks:ImproperInput','plotTracks: character input is not correct')
    %             end
    %         end
    %     end
end

%check inputs
if isempty(options.tracks)
    error('plotTracks:EmptyTracks','plotTracks: tracks field of options cannot not be empty')
end


% make figure if not indicated in input options
if isempty(options.h) && isempty(options.ha)
    h = figure; %make figure
else
    if ~isempty(options.ha)
        h = get(options.ha,'parent');
    else
        h = options.h; %get figure handle
    end
end

%set key press function
% if options.WindowKeyPressFcn
%     set(h,'WindowKeyPressFcn',@KeyPress);
% end

%setup axes
if isempty(options.ha)
    ha = findall(allchild(h),'type','axes'); %get axis handle
    if ~isempty(findall(ha,'tag','Colorbar'))
        ha = ha(ha ~= findall(ha,'tag','Colorbar'));
    end
    if isempty(ha)
        ha = axes('parent',h);
    end
    ha = ha(1);
else
    ha = options.ha;
end
options.ha = ha;
options.h = h;
if ~isempty(options.xlim) && ~isempty(options.ylim)
    set(ha,'plotboxaspectratio',[diff([options.xlim;options.ylim]')+1 max(diff([options.xlim;options.ylim]')+1)]) %set plotting aspect ratio
end

options = updateOptions(options);
[vert face] = setVertFace(options);

%make patch object
hp = makePatch(face,vert,'plotTrackPatch',options);

% create colorbar
if options.colorbar && ~isempty(options.minMaxVal) && options.minMaxVal(1)~=options.minMaxVal(2)
    colormap(options.h,options.cmap) %define colormap
    hcb = colorbar('peer',ha); %make colorbar
    set(get(hcb,'title'),'string',options.valueName) %set colorbar title
    set(findall(hcb,'tag','TMW_COLORBAR'),'ydata',options.minMaxVal) %set range for colorbar values
    set(hcb,'ylim',options.minMaxVal) %set range for colorbar figure
end
if ~isempty(options.view)
    view(options.ha,options.view)
end
resetAxes(options);

%setup data tip
dcm_obj = datacursormode(h);
set(dcm_obj,'updatefcn',@plotTracksDataCursor);
%add toolbar
% if options.WindowKeyPressFcn && isempty(findall(h,'userdata','plotTracksToolbar'))
%     addToolBar(h)
% end
% addToolBar(h)
% if isempty(findall(h,'tag','plotTracksMenu'))
%     addMenu(h)
% end
%add shadow on xy plane
hp = makeShadow(hp,options);


function resetAxes(options)
%reset data axes

%only change renderer and axis if not dipimage figure
drawnow
udataFig = get(options.h,'userdata');
if ~isfield(udataFig,'state')
    set(options.h,'renderer','zbuffer')
    axis(options.ha,'ij')
end
%     if ~isempty(options.view)
%         view(options.ha,options.view)
%     end
tightcheck = 0;
if ~isempty(options.xlim)
    xlim(options.ha,options.xlim);
    tightcheck = tightcheck+1;
end
if ~isempty(options.ylim)
    ylim(options.ha,options.ylim);
    tightcheck = tightcheck+1;
end
if isfield(options,'zlim') && ~isempty(options.zlim)
    zlim(options.ha,options.zlim);
    tightcheck = tightcheck+1;
end
if tightcheck ~= 3 && ~isfield(udataFig,'state')
    axis(options.ha,'tight')
end
xdiff = diff(get(options.ha,'xlim'));
ydiff = diff(get(options.ha,'ylim'));
zdiff = diff(get(options.ha,'zlim'));
if isfield(options,'DataAspectRatio')
    set(options.ha,'DataAspectRatio',options.DataAspectRatio)
else
    set(options.ha,'DataAspectRatio',[1 1 zdiff/max([xdiff ydiff])])
end
% drawnow

function [vert face] = setVertFace(options)
% find vertices with specified trackNum and in specified trange

if isfield(options,'zPlotFlag') && options.zPlotFlag
    vert = options.vert(:,[1 2 3]);
else
    vert = options.vert(:,[1 2 end-1]);
    vert(:,3) = (vert(:,3)-1)*options.t;
end
face = options.face(:,1:2);
vertIdx = ismember(options.vert(:,end),options.trackNum) & ismember(options.vert(:,end-1),options.trange);
% vert(:,3) = (vert(:,3)-1)*options.t;
vert(~vertIdx,:) = NaN;
face(any(~ismember(options.face(:,1:2),find(vertIdx)),2),:) = NaN;


function hp = makeShadow(hp,options)
%add shadow on xy plane

if nargin == 1
   options = get(hp,'userdata');
end

if options.addShadow && ~isempty(hp)
    shadOptions = options;
    [vert face] = setVertFace(shadOptions);
    shadOptions.addShadow = 0;
    % shadOptions.tagAppendix = 'Shadow';
    shadOptions.shadowParent = hp;
    %     if ~isempty(options.trange)
    %         shadOptions.trange = min(options.trange);
    %     else
    %         shadOptions.trange = 0;
    %     end
    shadOptions.colorByValue = 0;
    shadOptions.color = [.75 .75 .75];
    % shadOptions.face = options.face;
    shadOptions.markersize = options.markersize*options.addShadow;
    shadOptions.linewidth = options.linewidth*options.addShadow;
    %     shadOptions.vert(:,shadOptions.ndims+1) = shadOptions.trange;
    shadOptions = updateOptions(shadOptions);
    % vert(:,end) = shadOptions.trange;
    %     vert(:,end) = min(get(shadOptions.ha,'zlim'))-0.5;
    zdata = get(findall(shadOptions.h,'tag','plotTrackPatch'),'zdata');
    if iscell(zdata)
        zdata = cell2mat(zdata');
    end
    vert(:,end) = min(min(zdata))-shadOptions.t*0.5;
    
    %make patch object
    makePatch(face,vert,'plotTrackPatchShadow',shadOptions);
    % hp(2) = plotTracksV1(shadOptions);
    
    %repostion shadows z position
    allShadows = findall(shadOptions.h,'tag','plotTrackPatchShadow');
    faces = get(findall(shadOptions.h,'tag','plotTrackPatch'),'faces');
    zdata = get(findall(shadOptions.h,'tag','plotTrackPatch'),'zdata');
    if iscell(zdata)
        zdata = cell2mat(zdata')';
        faces = cell2mat(faces);
    else
        zdata = zdata';
    end
    for ii = 1:length(allShadows)
        udata = get(allShadows(ii),'userdata');
        [vert face] = setVertFace(udata);
        if any(any(~isnan(faces)))
            vert(:,end) = min(zdata(~isnan(faces)))-udata.t*1.5;
            makePatch(face,vert,'plotTrackPatchShadow',udata,allShadows(ii));
        end
    end
end


function options = updateOptions(options)
%update options for vertices, faces etc.

%get track numbers if not specified
if isempty(options.trackNum)
    options.trackNum = 1:size(options.tracks,1);
end
options.trackNum = options.trackNum(:)';
%get time range if not specified
if isempty(options.trange)
    options.trange = 1:size(options.tracks,2);
end
options.trange = options.trange(:)';

%compile vertices and faces for all tracks
if ~isfield(options,'updateOptions') || options.updateOptions
    options.ndims = size(options.tracks,3);
    %vertices for patch [positions t trackIdx]
    options.vert = zeros(sum(sum(logical(options.tracks(:,:,1)))),options.ndims+2);
    %faces for patch [vert1 vert2 trackIdx]
    options.face = zeros(size(options.vert,1)-size(options.tracks,1),3);
    %     %index for track that vertex belongs to
    %     options.vertTrackIdx = zeros(size(options.vert,1),1);
    %     %index for track that vertex belongs to
    %     options.vertFaceIdx = zeros(size(options.face,1),1);
    if ~isempty(options.tracksVal)
        options.val = zeros(size(options.vert,1),1);
    else
        options.val = [];
    end
    count = 1;
    count1 = 0;
    for ii = 1:size(options.tracks,1)
        nObsTrack = sum(logical(options.tracks(ii,:,1)));
        vertIdx = count+(0:nObsTrack-1);
        if nObsTrack
            if nObsTrack > 1
                options.vert(count+(0:nObsTrack-1),:) = [squeeze(options.tracks(ii,logical(options.tracks(ii,:,1)),:)) find(logical(options.tracks(ii,:,1)))' repmat(ii,[nObsTrack 1])];
            else
                options.vert(count+(0:nObsTrack-1),:) = [reshape(options.tracks(ii,logical(options.tracks(ii,:,1)),:),[1 options.ndims]) find(logical(options.tracks(ii,:,1)))' repmat(ii,[nObsTrack 1])];
                %         options.vert(count+(0:nObsTrack-1),end) = repmat(ii,[nObsTrack 1]);
            end
            options.face((count1+(1:nObsTrack-1)),:) = [vertIdx(1:end-1)' vertIdx(2:end)' repmat(ii,[nObsTrack-1 1])];
            %         options.vertFaceIdx(count+(1:nObsTrack-1)) = repmat(ii,[nObsTrack-1 1]);
            if ~isempty(options.tracksVal)
                options.val(count+(0:nObsTrack-1)) = squeeze(options.tracksVal(ii,logical(options.tracks(ii,:,1)))) ;
            end
            count = count+nObsTrack;
            count1 = count1+nObsTrack-1;
        end
    end
end

%get track coloration
if options.colorByValue
    if ~isfield(options,'updateOptions') || options.updateOptions
        if ~isempty(options.val) && any(options.val ~= 1)
            % get color information for plotting
            if isempty(options.minMaxVal)
                options.minMaxVal = [min(options.val) max(options.val)];
            end
            options.facevertexcdata = colorstretch(options.val,options.minMaxVal,options.cmap);
            if strcmp(options.linestyle,'-')
                options.edgecolor = 'interp';
            else
                options.edgecolor = 'flat';
            end
            options.markerfacecolor = 'flat';
            options.markeredgecolor = 'flat';
        else
            options.facevertexcdata = options.color;
            options.edgecolor = options.color;
            options.markerfacecolor = options.color;
            options.markeredgecolor = options.color;
        end
    end
else
    options.facevertexcdata = options.color;
    options.edgecolor = options.color;
    options.markerfacecolor = options.color;
    options.markeredgecolor = options.color;
end

options.updateOptions = 0;


function KeyPress(h,evnt)
% scroll feature

%variable inputs
if isempty(evnt)
    key = get(h,'userdata');
else
    key = evnt.Key;
end

% if strcmp(get(get(h,'parent'),'userdata'),'plotTracksToolbar') || strcmp(get(get(h,'parent'),'tag'),'plotTracksMenu')
%     hf = get(get(h,'parent'),'parent');
% else
%     hf = h;
% end

hf = ancestor(h,'figure');

% get userdata
switch key
    case 'r'
        resetTracks(hf);
    case 'h'
        syncSPTHSIplotFitResultsBrowser(hf);
    case 'd'
        initDataTip(hf);
    case 's'
        updatePlotTrackOptions(hf);
    case 'e'
        updateDipTrackOptions(hf);
end



function addMenu(h)
%add plotTrack menu to figure

hm = uimenu(h,'Label','plotTracks','Tag','plotTracksMenu');
uimenu(hm,'Label','reset tracks','Tag','resetMenu',...
    'Callback',@KeyPress,'userdata','r');
uimenu(hm,'Label','toggle datacursormode','Tag','dataTipMenu',...
    'Callback',@KeyPress,'userdata','d');
uimenu(hm,'Label','BrowseFitResults','Tag','SPTHSIplotFitResultsBrowserToggleMenu',...
    'Callback',@KeyPress,'userdata','h');
uimenu(hm,'Label','set plotTrackOptions','Tag','plotTracksOptionsMenu',...
    'Callback',@KeyPress,'userdata','s');
uimenu(hm,'Label','set dipTrackOptions','Tag','plotTracksOptionsMenu',...
    'Callback',@KeyPress,'userdata','e');



function addToolBar(h)
% add plotTracks toolbar to figure

% ht = uitoolbar(h);
set(h,'toolbar','figure')
ht = findall(h,'type','uitoolbar');
tags = get(get(ht,'children'),'tag');
% set(ht,'userdata','plotTracksToolbar')
if ~any(strcmp('SPTHSI:FitResultsBrowser',tags))
    load('cubeIcon','cdata');
    htt = uitoggletool(ht,'separator','on','cdata',cdata,'TooltipString','view box fitting in SPTHSI:FitResultsBrowser',...
        'clickedcallback',@KeyPress,'userdata','h','tag','SPTHSI:FitResultsBrowser');
end
if ~any(strcmp('resetPushTool',tags))
    load('colorlinesIcon','cdata');
    hpt = uipushtool(ht,'cdata',cdata,'TooltipString','reset tracks to original state',...
        'clickedcallback',@KeyPress,'userdata','r','tag','resetPushTool');
end
if ~any(strcmp('initDataTipPushTool',tags))
    load('colorlinesboxIcon','cdata');
    hpt = uipushtool(ht,'cdata',cdata,'TooltipString','initialize specialized data tip',...
        'clickedcallback',@KeyPress,'userdata','d','tag','initDataTipPushTool');
end

function hp = makePatch(face,vert,tag,options,hp)
%make Patch object using options

% if all(isnan(vert))
%     hp = [];
%     return;
% end

if isfield(options,'tagAppendix')
    tag = [tag options.tagAppendix];
end

if nargin == 5
    options = get(hp,'userdata');
    highlightTracksOriginal = options.highlightTracks;
    options.highlightTracks = 0;
    set(hp,'userdata',options);
    %     set(hp,'facevertexcdata', options.facevertexcdata, ...
    %         'edgecolor', options.edgecolor, 'facecolor', 'none','markersize',options.markersize,...
    %         'marker',options.marker,'markerfacecolor', options.markerfacecolor,...
    %         'linewidth',options.linewidth,'linestyle',options.linestyle,...
    %         'markeredgecolor',options.markeredgecolor);
    set(hp,'faces', face, 'vertices', vert, 'facevertexcdata', options.facevertexcdata, ...
        'edgecolor', options.edgecolor, 'facecolor', 'none','markersize',options.markersize,...
        'marker',options.marker,'markerfacecolor', options.markerfacecolor,...
        'linewidth',options.linewidth,'linestyle',options.linestyle,...
        'markeredgecolor',options.markeredgecolor);
     options.highlightTracks = highlightTracksOriginal;
     set(hp,'userdata',options);
else
    hp = patch('faces', face, 'vertices', vert, 'facevertexcdata', options.facevertexcdata, ...
        'edgecolor', options.edgecolor, 'facecolor', 'none','markersize',options.markersize,...
        'marker',options.marker,'markerfacecolor', options.markerfacecolor,...
        'linewidth',options.linewidth,'linestyle',options.linestyle,...
        'markeredgecolor',options.markeredgecolor,...
        'parent',options.ha,'tag',tag,'userdata',options);
end

if isfield(options,'displayName')
    set(hp,'displayname',options.displayName)
end


function hpNew = highlightTracks(hp,trackNum,linewidth,markersize,makeOthersTransparent,...
    linecolor,markeredgecolorgray,marker)
% highlight indicated track lines
%
% INPUTS
%   hp - patch handle
%   trackNum - vector with length N containing track numbers to highlight.
%   linewidth - linewidth for highlighted track. Default 4
%   markersize - marker size for highlighted track.
%   makeOthersTransparent - binary for graying out other trajectories. Default 1
%   linecolor - color for highligthed tracsk 1 by 3 vector for rgb
%   markeredgecolorgray - binary for making line marker edge color gray(1) or same color as
%                   track(0). Default 1

%check inputs
if ~exist('hp','var') || isempty(hp)
    return;
end
if ~exist('trackNum','var') ||isempty(trackNum)
    return;
end
if ~exist('linewidth','var') ||isempty(linewidth)
    linewidth = 4;
end
if ~exist('markersize','var') ||isempty(markersize)
    markersize = 6;
end
if ~exist('makeOthersTransparent','var') ||isempty(makeOthersTransparent)
    makeOthersTransparent = 1;
end
if ~exist('linecolor','var') ||isempty(linecolor)
    linecolor = [];
end
if ~exist('markeredgecolorgray','var') ||isempty(markeredgecolorgray)
    markeredgecolorgray = 1;
end
if ~exist('marker','var') ||isempty(marker)
    marker = 'o';
end

%find patch handle
if ~ishandle(hp)
    hptmp = findall(0,'tag','plotTrackPatch');
    hpIdx = round(hptmp - double(hp)) == 0;
    hp = hptmp(hpIdx);
end
if isempty(hp)
    return;
end

%find all existing highlighted tracks
hpHighlightOld = findall(get(hp,'parent'),'tag','plotTrackPatchHighlight');
trackNumOld = [];
if ~isempty(hpHighlightOld)
    for ii = 1:length(hpHighlightOld)
        udataOld = get(hpHighlightOld(ii),'userdata');
        if isfield(udataOld,'trackii')
            trackNumOld = [trackNumOld; udataOld.trackii];
        end
    end
end
delete(hpHighlightOld)
%get color information from figure
udata = get(hp,'userdata');
%find vertices for specifie
trackNum = unique([trackNumOld trackNum(:)']);

if makeOthersTransparent
    %change transparency of patch
    set(hp,'facevertexcdata', udata.facevertexcdata*.3+.7,...
        'linewidth',udata.linewidth*.75,'markersize',udata.markersize*.75);
end
%get marker edge color
if markeredgecolorgray
    if length(trackNum) > 1
        markeredgeC = repmat((.4:.4/(length(trackNum)-1):.8)',[1 3]);
    else
        markeredgeC = [.4 .4 .4];
    end
else if ~isempty(linecolor)
        markeredgeC = linecolor;
    else
        markeredgeC = [];
    end
end

%loop through tracks
if ~isempty(linecolor)
    udata.facevertexcdata = linecolor;
    udata.edgecolor = linecolor;
end
udata.markersize = markersize;
udata.linewidth = linewidth;
udata.marker = marker;
if size(markeredgeC,1) == length(trackNum)
    for ii = 1:length(trackNum)
        udata.trackii = trackNum(ii);
        % find vertices with specified trackNum and in specified trange
        %         vert = udata.vert(:,[1 2 end-1]);
        %         face = udata.face(:,1:2);
        %         vertIdx = udata.vert(:,end) == udata.trackii & ismember(udata.vert(:,end-1),udata.trange);
        %         face(any(~ismember(udata.face(:,1:2),find(vertIdx)),2),:) = [];
        [vert face] = setVertFace(udata);
        vertIdx = udata.vert(:,end) == udata.trackii & ismember(udata.vert(:,end-1),udata.trange);
        face(any(~ismember(udata.face(:,1:2),find(vertIdx)),2),:) = [];
        %set marker edge color for track ii
        udata.markeredgecolor = markeredgeC(ii,:);
        %make patch
        hpNew = makePatch(face,vert,'plotTrackPatchHighlight',udata);
    end
else
    udata.trackNum = trackNum;
    [vert face] = setVertFace(udata);
    hpNew = makePatch(face,vert,'plotTrackPatchHighlight',udata);
end

function initDataTip(h)
% initialize specialized data tip
%
% h - figure handle

dcm_obj = datacursormode(h);
if isempty(get(dcm_obj,'updatefcn'))
    set(dcm_obj,'updatefcn',@plotTracksDataCursor);
    disp('plotTracks: specialized data tip initialized')
end
if strcmp(dcm_obj.Enable,'on')
    datacursormode off
    disp('plotTracks: datacursormode toggled off')
else
    datacursormode on
    disp('plotTracks: datacursormode toggled on')
end



function hp = resetTracks(h)
% reset patch
%
% INPUT
%   h - figure handle
% OUTPUT
%   hp - updated patch handles

%delete all plotTrackPatchHighlight
hpH = findall(h,'type','patch','tag','plotTrackPatchHighlight');
delete(hpH)
%find all plotTrackPatchShadow
hpS = findall(h,'type','patch','tag','plotTrackPatchShadow');
shadowParent = [];
for ii = 1:length(hpS)
    udataShadow = get(hpS(ii),'userdata');
    shadowParent(ii) = udataShadow.shadowParent;
end
%find all plotTrackPatch
hp = findall(h,'type','patch','tag','plotTrackPatch');
% if length(hp)>1
%     idx = find(~cellfun('isempty',strfind(get(hp,'tag'),'plotTrackPatch')))';
% else
%     idx = ~isempty(strfind(get(hp,'tag'),'plotTrackPatch'));
% end
for ii = 1:length(hp)
    udata = get(hp(ii),'userdata');
    delete(hp(ii))
    udata.updateOptions = 1;
    udata = updateOptions(udata);
    [vert face] = setVertFace(udata);
    shadowIdx = hp(ii) == shadowParent;
    hp(ii) = makePatch(face,vert,'plotTrackPatch',udata);
    %make new shadow
    if any(shadowIdx)
        delete(hpS(shadowIdx))
        resetAxes(udata);
    end
    makeShadow(hp(ii));
end

%update dipTrack
udata = get(h,'userdata');
if ~isempty(udata) && isfield(udata,'curslice')
    if udata.curslice == udata.maxslice
        dipmapping(h,'slice',udata.curslice+-1);
    else
        dipmapping(h,'slice',udata.curslice+1);
    end
    dipmapping(h,'slice',udata.curslice);
end
hp = findall(h,'type','patch','tag','plotTrackPatch');





function updatePlotTrackOptions(h)
%update plotTracks options

hp = findall(h,'tag','plotTrackPatch');
updateFields = {'colorByValue' 'linewidth' 'linestyle' 'valueName'...
    'trackNum' 'view' 'color' 'minMaxVal' 'markersize'...
    'marker' 'trange' 'addShadow' 'highlightTracks' 'displayName' 'cmap'};
%reset
% resetTracks(h);
for ii = 1:length(hp)
    udata(ii) = get(hp(ii),'userdata');
end

for ii = 1:length(hp)
    %     udata = get(hp(ii),'userdata');
    %highlight patch
    highlightTracks(hp(ii),1:size(udata(ii).tracks,1),[],[],[],[],0);
    %Set fields with gui
    for jj = 1:length(updateFields)
        udataTmp.(updateFields{jj}) = udata(ii).(updateFields{jj});
    end
    udataTmp = StructDlg(udataTmp,sprintf('Set plotTrackOptions for group %i of %i',ii, length(hp)));
    if ~isempty(udataTmp)
        for jj = 1:length(updateFields)
            udata(ii).(updateFields{jj}) = udataTmp.(updateFields{jj});
        end
    end
    udata(ii) = updateOptions(udata(ii));
    set(hp(ii),'userdata',udata(ii))
    %reset
    hp = resetTracks(h);
    %     hp = findall(h,'tag','plotTrackPatch');
end


function updateDipTrackOptions(h)
%update plotTracks options

udata = get(h,'userdata');
if ~isempty(udata)
    updateFields = {'headMarker' 'headMarkersize' 'npoints'};
    %set headMarker field if not defined (back compatibility)
    if ~isfield(udata.dipTrackData,'headMarker') || isempty(udata.dipTrackData.headMarker)
        udata.dipTrackData.headMarker = 'o';
    end
    %Set fields with gui
    for jj = 1:length(updateFields)
        dipTrackDataTmp.(updateFields{jj}) = udata.dipTrackData.(updateFields{jj});
    end
    dipTrackDataTmp = StructDlg(dipTrackDataTmp,'Set dipTrackOptions');
    if ~isempty(dipTrackDataTmp)
        for jj = 1:length(updateFields)
            udata.dipTrackData.(updateFields{jj}) = dipTrackDataTmp.(updateFields{jj});
        end
    end
    set(h,'userdata',udata);
    resetTracks(h);
    %     if udata.curslice == udata.maxslice
    %         dipmapping(h,'slice',udata.curslice+-1);
    %     else
    %         dipmapping(h,'slice',udata.curslice+1);
    %     end
    %     dipmapping(h,'slice',udata.curslice);
end



function syncSPTHSIplotFitResultsBrowser(h)
% toggle syncing SPTHSIplotFitResultsBrowser with data tip
%
% h - figure handle

%get toggle tool handle
htt = findall(h,'tag','SPTHSI:FitResultsBrowser');%check for SPTHSI_FitResultsBrowser and data linkage
if isempty(htt)
    addToolBar(h)
    htt = findall(h,'tag','SPTHSI:FitResultsBrowser');%check for SPTHSI_FitResultsBrowser and data linkage
end
set(h,'Interruptible','off','BusyAction','cancel')
hp = findall(h,'tag','plotTrackPatch');
for ii = 1:length(hp)
    if get(htt,'state') %get toggle tool state
        udata = get(hp(ii),'userdata');
        try
            %             if ~isfield(udata,'sptBrowserObj')
            %                 disp('plotTrack: loading SPT object..');
            %                 tic
            %                 sptBrowserObj = SPTHSI_FitResultsBrowser(SPTHSI.loadFile);
            %                 if ii > 1
            %                     assignin('base',sprintf('sptBrowserObj%i',ii),sptBrowserObj);
            %                     fprintf('sptBrowserObj%i variable created in workspace\n',ii)
            %                 else
            %                     assignin('base','sptBrowserObj',sptBrowserObj);
            %                     fprintf('sptBrowserObj variable created in workspace')
            %                 end
            %                 udata.sptBrowserObj = sptBrowserObj;
            %                 %                 udata.sptBrowserObj = SPTHSI_FitResultsBrowser(sptBrowserObj);
            %                 toc
            %             end
            udata.sptBrowserObj.SptObj.checkDataLink;
            udata.sptBrowserObj.BrowseFlag = 1;
            set(hp(ii),'userdata',udata);
        catch me
            try
                %look for spt file 2 folders up
                [pathName, baseFigName] = fileparts(get(h,'FileName'));
                tmpIdx = strfind(baseFigName,'-');
                baseName = baseFigName(1:tmpIdx(end)-1);
                sptFile = fullfile(fileparts(fileparts(pathName)),[baseName '.spt']);
                udata.sptBrowserObj.SptObj = SPT.loadFile(sptFile);
                udata.sptBrowserObj.SptObj.checkDataLink;
                udata.sptBrowserObj.BrowseFlag = 1;
                set(hp(ii),'userdata',udata);
            catch me
                try
                    %look for spt file 3 folders up
                    sptFile = fullfile(fileparts(fileparts(fileparts(pathName))),[baseName '.spt']);
                    udata.sptBrowserObj.SptObj = SPT.loadFile(sptFile);
                    udata.sptBrowserObj.SptObj.checkDataLink;
                    udata.sptBrowserObj.BrowseFlag = 1;
                    set(hp(ii),'userdata',udata);
                catch me
                    if isfield(udata,'sptBrowserObj')
                        udata = rmfield(udata,'sptBrowserObj');
                        set(hp(ii),'userdata',udata);
                    end
                    set(htt,'state','off') %set toggle tool state
                    warning('plotTrack:NoSptPlotObj','plotTrack: synchronization with sptPlot.FitResultsBrowser not available for current plot');
                end
            end
            
        end
    else
        if isfield(udata,'sptBrowserObj')
            udata.sptBrowserObj.BrowseFlag = 1;
        end
    end
    
end


function output_txt = plotTracksDataCursor(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).


pos = get(event_obj,'Position');
output_txt = {['X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)]};

try
    % If there is a Z-coordinate in the position, display it as well
    if length(pos) > 2
        output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
    end
    
    %get patch handle
    hp = get(event_obj,'target');
    if ~strcmp(get(hp,'type'),'patch') || (isempty(strfind(get(hp,'tag'),'plotTrackPatch'))...
            && isempty(strfind(get(hp,'tag'),'dipTrackLocMarker')))
        output_txt = [];
        return;
    end
    
    displayName = get(hp,'displayname');
    if ~isempty(displayName)
        output_txt{end+1} = sprintf('%s', displayName);
    end
    
    udata =  get(hp,'userdata');
    % find vertex index for current point
    idx = get(event_obj,'dataindex');
    % display trackNum
    trackNum = udata.vert(idx,end);
    output_txt{end+1} = sprintf('trackNum: %i', trackNum);
    % If there is a Val field in userdata, display it as well
    if isfield(udata,'val') && ~isempty(udata.val)
        output_txt{end+1} = sprintf('Val: %.3g', udata.val(idx));
        output_txt{end+1} = sprintf('meanVal: %.3g', mean(udata.val(udata.vert(:,end) == udata.vert(idx,end))));
    end
    
    % If there is a Val field in userdata, display it as well
    zIdx = udata.vert(idx,end-1)/udata.t;
    if isfield(udata,'subRegionModel') && ~isempty(udata.subRegionModel)
        output_txt{end+1} = sprintf('boxIdxA: %i', udata.subRegionModel(trackNum,zIdx,1));
        output_txt{end+1} = sprintf('boxIdxF: %i', udata.subRegionModel(trackNum,zIdx,2));
        output_txt{end+1} = sprintf('model: %i', udata.subRegionModel(trackNum,zIdx,3));
        output_txt{end+1} = sprintf('idx: %i', udata.subRegionModel(trackNum,zIdx,4));
    end
    
    % If there is a ID field in userdata, display it
    if isfield(udata,'ID') && ~isempty(udata.ID)
        if ischar(udata.ID)
            output_txt{end+1} = sprintf('ID: %s', udata.ID);
        else if isnumerical(udata.ID)
                output_txt{end+1} = sprintf('ID: %g', udata.ID);
            end
        end
    end
    
    %highlight trackLines
    % highlightTracks(get(event_obj,'target'),udata.vert(idx,end))
    
    %highlight track
    h = get(get(get(event_obj,'target'),'parent'),'parent');
    dcm_obj = datacursormode(h);
    info_cell = squeeze(struct2cell(getCursorInfo(dcm_obj)));
    hLines = findall(cell2mat(info_cell(1,:)),'tag','plotTrackPatchHighlight');
    if isempty(dcm_obj)
        resetTracks(h);
    else
        if udata.highlightTracks
            hpHighlight = findall(hp,'tag','plotTrackPatchHighlight');
            trackNum = udata.vert(idx,end);
            for ii = 1:length(hpHighlight)
                udataHighlight = get(hpHighlight(ii),'userdata');
                trackNum = [trackNum;udataHighlight.trackii];
            end
            try
                highlightTracks(hp,trackNum);
            catch me
                output_txt{end+1} = sprintf('highlightTracks error\nmessage:%s\n',...
                    me.stack(end).file,me.stack(end).name,me.stack(end).line,me.message);
            end
        end
    end
    if isfield(udata,'sptBrowserObj') && isfield(udata,'subRegionModel')
        %get toggle tool handle
        htt = findall(h,'tag','SPTHSI:FitResultsBrowser');%check for sptBrowserObj and data linkage
        if isa(udata.sptBrowserObj,'SPTHSI_FitResultsBrowser') ...
                && udata.sptBrowserObj.BrowseFlag
            %                 && strcmp(get(htt,'state'),'on') ...
            %                 && strcmp(get(dcm_obj,'Enable'),'on')
            if zIdx ~= udata.sptBrowserObj.FrameIdx ...
                    && udata.subRegionModel(trackNum,zIdx,2) == udata.sptBrowserObj.BoxIdx
                udata.sptBrowserObj.FrameIdx = zIdx;
            else if zIdx == udata.sptBrowserObj.FrameIdx ...
                    && udata.subRegionModel(trackNum,zIdx,2) ~= udata.sptBrowserObj.BoxIdx
                udata.sptBrowserObj.BoxIdx = udata.subRegionModel(trackNum,zIdx,2);
                else if zIdx ~= udata.sptBrowserObj.FrameIdx ...
                    && udata.subRegionModel(trackNum,zIdx,2) ~= udata.sptBrowserObj.BoxIdx
                udata.sptBrowserObj.BrowseFlag = 0;
                udata.sptBrowserObj.FrameIdx = zIdx;
                udata.sptBrowserObj.BrowseFlag = 1;
                udata.sptBrowserObj.BoxIdx = udata.subRegionModel(trackNum,zIdx,2);
                %                 %                 disp(['plotTracks: updating frame ' num2str(zIdx) ' box ' num2str(udata.subRegionModel(trackNum,zIdx,2)) ' for SPTHSIplotFitResultsBrowser...'])
                %                 %                 tic
                %                 udata.sptBrowserObj.SptObj.plotFitResultsBrowser(zIdx,udata.subRegionModel(trackNum,zIdx,2));
                %                 %                 toc
                    end
                end
            end
        end
    end
    
catch me
    output_txt{end+1} = sprintf('plotTracksDataCursor error\nmessage:%s',...
        me.stack(end).file,me.stack(end).name,me.stack(end).line,me.message);
end

