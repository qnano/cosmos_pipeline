function varargout = dipTrack(varargin)
% dipTrack     dipimage tracking figure
%
% USAGE:
%   h = dipTrack(options); %setup dipTrack figure
%   h = dipTrack(fh); %setup listener for dipTrack figure
%



if nargin > 2
    error('dipTrack:ToManyInputs','dipTrack: 0 to 2 inputs required');
end

switch nargin
    case 0
        options = dipTrackSetOptions;
    case 1
        if ishandle(varargin{1})
            h = varargin{1};
            udata = get(h,'userdata');
            if isfield(udata,'dipTrackData')
                %add back compatibility for dipTrack figures made with plotTracks
                if ~isfield(udata.dipTrackData,'plotTracksOptions')
                    udata.dipTrackData.plotTracksOptions = plotTracksV1SetOptions;
                    udata.dipTrackData.plotTracksOptions.marker = '.';
                    udata.dipTrackData.plotTracksOptions.tracks = udata.dipTrackData.tracks;
                    if size(udata.dipTrackData.c,1) == size(udata.dipTrackData.plotTracksOptions.tracks,1)
                        udata.dipTrackData.plotTracksOptions.tracksVal = repmat((1:size(udata.dipTrackData.c,1))',[1 size(udata.dipTrackData.plotTracksOptions.tracks,2)]);
                        udata.dipTrackData.plotTracksOptions.minMaxVal = [min(udata.dipTrackData.plotTracksOptions.tracksVal(:)) max(udata.dipTrackData.plotTracksOptions.tracksVal(:))];
                        udata.dipTrackData.plotTracksOptions.cmap = udata.dipTrackData.c;
                    end
                    if iscell(udata.dipTrackData.trackNum)
                        trackNum = cell2mat(udata.dipTrackData.trackNum);
                    else
                        trackNum = udata.dipTrackData.trackNum;
                    end
                    udata.dipTrackData.plotTracksOptions.colorbar = 0;
                    udata.dipTrackData.plotTracksOptions.WindowKeyPressFcn = 0;
                    udata.dipTrackData.plotTracksOptions.view = [];
                    %                     udata.dipTrackData.plotTracksOptions.addShadow = 0;
                    udata.dipTrackData.plotTracksOptions.trackNum = trackNum;
                    udata.dipTrackData.plotTracksOptions.linewidth = udata.dipTrackData.linewidth(2);
                    udata.dipTrackData.plotTracksOptions.h = h;
                    udata.dipTrackData.plotTracksOptions.ha = findall(h,'type','axes');
                    udata.dipTrackData.headMarkersize = udata.dipTrackData.markersize;
                    udata.dipTrackData.trackToggle = 1;
                    %delete all prior tracks and localization markers
                    delete(findall(udata.dipTrackData.plotTracksOptions.ha,'tag','trackLine'));
                    delete(findall(udata.dipTrackData.plotTracksOptions.ha,'tag','dipTrackLocMarker'));
                end
                udata.dipTrackData.h = h;
                addlistener(udata.dipTrackData.h,'UserData','PostSet',@curslice);
                set(h,'userdata',udata);
                set(h,'KeyPressFcn',@dipTrack);
                return;
            else
                error('dipTrack:InputMustBeDipTrack','dipTrack: input figure must be initialized using dipTrack')
            end
        end
        options = varargin{1};
    case 2
        %dipimage keypress call back
        dipshow('DIP_callback','KeyPressFcn')
        h = varargin{1};
        udata = get(h,'UserData');

        switch varargin{2}.Key
            case {'t','T'} % initialize trackTable
                %                 update_trackTable(h,udata);
            case {'r','R'} % toggle track plotting
                udata.dipTrackData.trackToggle = ~udata.dipTrackData.trackToggle;
                update_tracks(udata);
                set(h,'userdata',udata);
            case {'d','D'}
                plotTracksV1('initDataTip(gcf)');
        end
        return
end
%initialize dipData
dipTrackData.headMarkersize = options.headMarkersize;
dipTrackData.dipTrackObj = dipTrackObj(0);
dipTrackData.trackTableToggle = 0;
dipTrackData.trackToggle = 1;
dipTrackData.npoints = options.npoints;
%initilize figure
if isempty(options.h)
    h = figure;
else
    h = options.h;
end
dipTrackData.h = h;
%show image
if options.updateImage
    if isempty(options.im)
        warning('dipTrack:NoIm','dipTrack: options.im is empty initializing with newim(256,256,10)')
        if ~isempty(options.plotTracksOptions.tracks)
            sz=round(max(dip_image(options.plotTracksOptions.tracks)));
            dipshow(h,newim([sz,sz,size(options.plotTracksOptions.tracks,2)]));
        else
            dipshow(h,newim([256,256,10]));
        end
    else
        dipshow(h,options.im);
    end
end
%set values specific to plotTracksOptions
for ii = 1:length(options.plotTracksOptions)
    dipTrackData.plotTracksOptions(ii) = options.plotTracksOptions(ii);
    dipTrackData.plotTracksOptions(ii).h = h;
    dipTrackData.plotTracksOptions(ii).ha = findall(h,'type','axes');
    dipTrackData.plotTracksOptions(ii).colorbar = 0;
    %     dipTrackData.plotTracksOptions(ii).WindowKeyPressFcn = 0;
    dipTrackData.plotTracksOptions(ii).view = [];
    %     dipTrackData.plotTracksOptions(ii).addShadow = 0;
    if strcmp(dipTrackData.plotTracksOptions(ii).marker,'none')
        dipTrackData.plotTracksOptions(ii).marker = '.';
    end
end
%add listener
dipTrackData.lh = addlistener(h,'UserData','PostSet',@curslice);
%add dipTrackData to userdata
udata = get(h,'userdata');
udata.dipTrackData = dipTrackData;
set(h,'userdata',udata);
%change KeyPressFcn
set(h,'KeyPressFcn',@dipTrack,'busyaction','cancel','interruptible','off');
%add listener
% addlistener(udata.dipTrackData.h,'UserData','PostSet',@curslice);

%initialize tracks
dipmapping(h,'slice',1);
dipmapping(h,'slice',0);

% %add initialization menu
% if isempty(findall(h,'Tag','initializeMenu'))
%     hm = findall(h,'Label','plotTracks','Tag','plotTracksMenu');
%     uimenu(hm,'Label','initialize dipTrack','Tag','initializeMenu',...
%         'Callback',@dipTrackInitialize);
% end
% ht = findall(h,'type','uitoolbar');
% tags = get(get(ht,'children'),'tag');
% if ~any(strcmp('initDipTrackPushTool',tags))
%     load(fullfile(matlabroot,'/toolbox/matlab/icons/','arrow.mat'),'arrowCData');
%     hpt = uipushtool(ht,'cdata',arrowCData,'TooltipString','initialize dipTrack',...
%         'clickedcallback',@KeyPress,'userdata','d','tag','initDipTrackPushTool');
% end


%output figure handle
if nargout == 1
    varargout{1} = h;
end


function dipTrackInitialize(h,~)
%initialize dipTrack

hf = get(get(h,'parent'),'parent');
dipTrack(hf);


function update_tracks(udata)
%update trajectories in dipTrack figure

if isempty(udata.slicing)
    return;
end

% plotTracksOptions = udata.dipTrackData.plotTracksOptions;
hp = findall(udata.dipTrackData(1).plotTracksOptions(1).ha,'tag','plotTrackPatch');
if isempty(hp)
    plotTracksOptions = udata.dipTrackData.plotTracksOptions;
else
    for ii = 1:length(hp)
        plotTracksOptions(ii) = get(hp(ii),'userdata');
    end
end

%delete all prior tracks and localization markers
delete(findall(udata.dipTrackData(1).plotTracksOptions(1).ha,'tag','plotTrackPatch'));
delete(findall(udata.dipTrackData(1).plotTracksOptions(1).ha,'tag','plotTrackPatchShadow'));
delete(findall(plotTracksOptions(1).ha,'tag','plotTrackPatchHighlight'));
delete(findall(plotTracksOptions(1).ha,'tag','dipTrackLocMarker'));

if udata.dipTrackData.trackToggle
    if udata.dipTrackData.npoints > udata.maxslice+1
        npoints = udata.maxslice+1;
    else
        npoints = udata.dipTrackData.npoints;
    end
    trange = udata.curslice+(1:-1:-npoints);
    trange(trange<0) = [];
    for ii = 1:length(plotTracksOptions)
        plotTracksOptions(ii).trange = trange;
        [h ha hp] = plotTracksV1(plotTracksOptions(ii));
        if ~isempty(hp)
            plotTracksOptions1 = get(hp,'userdata');
            plotLocMarker(udata,plotTracksOptions1);
        end
    end
end


function plotLocMarker(udata,options)
%plot localization marker

tag = 'dipTrackLocMarker';
if isfield(udata.dipTrackData,'headMarker') && length(udata.dipTrackData.headMarker) == length(options.trackNum)
    for ii = 1:length(options.trackNum)
        if options.tracks(options.trackNum(ii),udata.curslice+1,1)
            line(options.tracks(options.trackNum(ii),udata.curslice+1,1),...
                options.tracks(options.trackNum(ii),udata.curslice+1,2),...
                udata.curslice,'tag',tag,'marker',udata.dipTrackData.headMarker(ii),...
                'markersize',udata.dipTrackData.headMarkersize(1),...
                'markerfacecolor',[1 1 1],...
                'markeredgecolor',[0 0 0])
            if isfield(udata.dipTrackData,'headColor')
                line(options.tracks(options.trackNum(ii),udata.curslice+1,1),...
                    options.tracks(options.trackNum(ii),udata.curslice+1,2),...
                    udata.curslice,'tag',tag,'marker',udata.dipTrackData.headMarker(ii),...
                    'markersize',udata.dipTrackData.headMarkersize(2),...
                    'markerfacecolor',udata.dipTrackData.headColor(ii,:),...
                    'markeredgecolor',[1 1 1])
            else
                line(options.tracks(options.trackNum(ii),udata.curslice+1,1),...
                    options.tracks(options.trackNum(ii),udata.curslice+1,2),...
                    udata.curslice,'tag',tag,'marker',udata.dipTrackData.headMarker(ii),...
                    'markersize',udata.dipTrackData.headMarkersize(2),...
                    'markerfacecolor',colorstretch(options.tracksVal(options.trackNum(ii),udata.curslice+1,1),options.minMaxVal,options.cmap),...
                    'markeredgecolor',[1 1 1])
            end
        end
    end
else
    if isfield(udata.dipTrackData,'headMarker') && ~isempty(udata.dipTrackData.headMarker)
        options.marker = udata.dipTrackData.headMarker(1);
    else
        options.marker = 'o';
    end
    facevertexcdata = options.facevertexcdata;
    %     options.marker = 'o';
    options.markersize = udata.dipTrackData.headMarkersize(1);
    options.facevertexcdata = facevertexcdata*0.5;
    options.markeredgecolor = [1 1 1];
    options.trange = udata.curslice+1;
    vert = options.vert(:,[1 2 end-1]);
    face = options.face(:,[1 2]);
    
    vertIdx = ismember(options.vert(:,end),options.trackNum) & ismember(options.vert(:,end-1),options.trange);
    vert(~vertIdx,:) = NaN;
    face(all(~ismember(face(:,1:2),find(vertIdx)),2),:) = NaN;
    %     face(~ismember(face(:,2),find(vertIdx)),:) = NaN;
    if options.trange-1 == udata.maxslice
        face = face(:,[2 2]);
    else
        face = face(:,[1 1]);
    end
    %     if options.trange == 1
    %         face = face(:,[1 1]);
    %     else
    %         face = face(:,[2 2]);
    %     end

    hp = patch('faces', face, 'vertices', vert, 'facevertexcdata', options.facevertexcdata, ...
        'edgecolor', options.edgecolor, 'facecolor', 'none','markersize',options.markersize,...
        'marker',options.marker,'markerfacecolor', options.markerfacecolor,...
        'linewidth',options.linewidth,'linestyle',options.linestyle,...
        'markeredgecolor',options.markeredgecolor,...
        'parent',options.ha,'tag',tag,'userdata',options);
    
    options.markersize = udata.dipTrackData.headMarkersize(2);
    options.facevertexcdata = facevertexcdata;
    if size(facevertexcdata,1) == 1
        options.markeredgecolor = facevertexcdata;
    else
        options.markeredgecolor = 'flat';
    end
    hp = patch('faces', face, 'vertices', vert, 'facevertexcdata', options.facevertexcdata, ...
        'edgecolor', options.edgecolor, 'facecolor', 'none','markersize',options.markersize,...
        'marker',options.marker,'markerfacecolor', options.markerfacecolor,...
        'linewidth',options.linewidth,'linestyle',options.linestyle,...
        'markeredgecolor',options.markeredgecolor,...
        'parent',options.ha,'tag',tag,'userdata',options);
end

function curslice(h,evnt)
%check updating of curslice field in userdata 

%get userdata
udata = get(evnt,'NewValue');
%get current slice
if isfield(udata,'curslice') && isfield(udata,'dipTrackData') &&... 
        udata.dipTrackData.dipTrackObj.slice ~= udata.curslice
    udata.dipTrackData.dipTrackObj.slice = udata.curslice; %update slice property
    update_tracks(udata); %update figure
    %     set(udata.dipTrackData.h,'userdata',udata)
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


function update_trackTable(h,udata)
%update trackTable

%populate fields for backcompatibility
if ~isfield(udata.dipTrackData,'ntrackTable')
    udata.dipTrackData.ntrackTable = 5;
end
if ~isfield(udata.dipTrackData,'trackTableDeleteLast')
    udata.dipTrackData.trackTableDeleteLast = 1;
end
if ~isfield(udata.dipTrackData,'trackTableShowSummary')
    udata.dipTrackData.trackTableShowSummary = 1;
end

locs=dipgetcoords(h,1);

if udata.dipTrackData.trackTableDeleteLast
    curFigs = get(findall(0,'tag','dipTrackTable'),'parent');
    if iscell(curFigs)
        curFigs = cell2mat(curFigs);
    end
    delete(curFigs)
end

%get track numbers
if iscell(udata.dipTrackData.trackNum)
    trackNum = [];
    for ii = 1:numel(udata.dipTrackData.trackNum)
        trackNum = [trackNum;udata.dipTrackData.trackNum{ii}(:)];
    end
else
    if isempty(udata.dipTrackData.trackNum)
        trackNum = 1:size(udata.dipTrackData.tracks,1);
    end
end
%get tracks observed in the current frame
trackNum = trackNum(logical(udata.dipTrackData.tracks(trackNum,locs(3)+1,1)));
%sort tracks by distance from locs
dist = sum((repmat(reshape(locs(1:2),[1 1 2]),[length(trackNum) 1 1])-udata.dipTrackData.tracks(trackNum,locs(3)+1,:)).^2,3);
[val idx] = sort(dist);
trackNum = trackNum(idx);
%get track information
trackNum = trackNum(1:min(udata.dipTrackData.ntrackTable,length(trackNum)));
x = udata.dipTrackData.tracks(trackNum,locs(3)+1,1);
y = udata.dipTrackData.tracks(trackNum,locs(3)+1,2);
c = round(udata.dipTrackData.c(trackNum,:)*255);
% c = zeros(length(trackNum),3);
% obs = zeros(length(trackNum),3); %first observations,last observation, number of total observations
data = cell(length(trackNum),10);
for ii = 1:length(trackNum)
    obsfirst = find(udata.dipTrackData.tracks(trackNum(ii),:,1),1,'first');
    obslast = find(udata.dipTrackData.tracks(trackNum(ii),:,1),1,'last');
    obsall = sum(logical(udata.dipTrackData.tracks(trackNum(ii),:,1)));
    %     data(ii,:) = {trackNum(ii) x(ii) y(ii) obsfirst obslast obsall false sprintf('%0.1g ',c(ii,:)) '' ''};
    cii = sprintf('%s%i%s%i%s%i%s%i%s%i%s%i%s','<html><span style="background-color: rgb(',...
        c(ii,1),',',c(ii,2),',',c(ii,3),');">(',...
        c(ii,1),',',c(ii,2),',',c(ii,3),')</span></html>');
    data(ii,:) = {trackNum(ii) x(ii) y(ii) obsfirst obslast obsall false cii '' ''};
end
cnames = {'track#' 'x' 'y' 'First obs' 'Last obs' '# obs','append' 'color' 'identifier' 'comments'};
cformat= {'numeric' 'numeric' 'numeric' 'numeric' 'numeric' 'numeric' 'logical' 'char' 'char' 'char'};
ceditable = logical([0 0 0 0 0 0 1 0 1 1]);
th = figure;
namestr = sprintf('%s%i%s%i%s%i','dipTrackTable(x-',locs(1),',y-',locs(2),',t-',locs(3),')');
t=uitable(th,'ColumnName', cnames,'ColumnFormat', cformat, 'Data', data,...
    'columneditable', ceditable,'tag','dipTrackTable');
%update table
pos=get(t,'extent');
set(t, 'Position', pos );
udata1.dipTrackData.trackTableShowSummary = udata.dipTrackData.trackTableShowSummary;
set(th,'Position',pos+[10 75 0 0],'name',namestr,...
    'DeleteFcn',@deleteTrackTable,'userdata',udata1);
set(h,'userdata',udata)
figure(h);

%insert code for closing fcn
%add export fcn. export trackTable to workspace and append

function deleteTrackTable(fh,evnt)


%get table data
dipTrackTableData = get(findall(fh,'tag','dipTrackTable'),'data');
%find data to append
if ~isempty(dipTrackTableData)
    appendIdx = cell2mat(dipTrackTableData(:,7));
    dipTrackTableData = dipTrackTableData(appendIdx,:);
    %eliminate append column
    dipTrackTableData(:,7) = [];
end
%update summary table
update_trackTableSummary(fh,dipTrackTableData)



function update_trackTableSummary(fh,dipTrackTableData)

%get current figure handle
curSummaryTable = findall(0,'tag','diptracktablesummary');
if ~isempty(curSummaryTable)
    dipTrackTableData = [get(curSummaryTable,'data');dipTrackTableData];
else if evalin('base','exist(''dipTrackTableData'',''var'')')
        dipTrackTableData = [evalin('base','dipTrackTableData');dipTrackTableData];
    end
end
udata = get(fh,'userdata');
%make summary table
curfigs = get(curSummaryTable,'parent');
if iscell(curfigs)
    curfigs = cell2mat(curfigs);
end
delete(curfigs)
if ~isempty(dipTrackTableData)
    if udata.dipTrackData.trackTableShowSummary
        cnames = {'track#' 'x' 'y' 'first obs' 'last obs' '# obs' 'color' 'identifier' 'comments' 'delete'};
        cformat= {'numeric' 'numeric' 'numeric' 'numeric' 'numeric' 'numeric' 'char' 'char' 'char' 'logical'};
        ceditable = logical([0 0 0 0 0 0 0 1 1 1]);
        th = figure;
        namestr = sprintf('%s%i%s%i%s%i','dipTrackTable(summary)');
        t=uitable(th,'columnname', cnames,'columnformat', cformat, 'data', dipTrackTableData,...
            'columneditable', ceditable,'tag','diptracktablesummary',...
            'celleditcallback',@update_dipTrackTableSummary);
        monpos = get(0,'monitorposition');
        [val idx] = min(sum(monpos(:,1:2),2)); %find lower left monitor
        pos=get(t,'extent');
        if pos(4) > monpos(idx,4)-100
            pos(4) = monpos(idx,4)-150;
        end
        if pos(3) > monpos(idx,3)-100
            pos(3) = monpos(idx,3)-150;
        end
        pos1 = [monpos(idx,[3 4]) 0 0]+[-pos([3 4])-75 pos([3 4])];
        set(t, 'position', pos );
        set(th,'position',pos1,'userdata',udata,'name',namestr);
        
    end
end
assignin('base','dipTrackTableData',dipTrackTableData)


function update_dipTrackTableSummary(t,evnt)

if evnt.Indices(2) == 10
    button = questdlg(sprintf('%s%i%s','Do you want to delete row',evnt.Indices(1),'?'),...
        'dipTrackTable(summary)','YES','NO','YES');
    switch button
        case 'YES'
            th = get(t,'parent');
            %             dipTrackTableData = evalin('base','dipTrackTableData');
            curSummaryTable = findall(0,'tag','diptracktablesummary');
            dipTrackTableData = get(curSummaryTable,'data');
            dipTrackTableData(evnt.Indices(1),:) = [];
            assignin('base','dipTrackTableData',dipTrackTableData);
            set(curSummaryTable,'data',dipTrackTableData);
            %             update_trackTableSummary(th,[]);
            return
        case 'NO'
            return
    end
end

