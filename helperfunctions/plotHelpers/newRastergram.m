function [evi,cev,cia,meanTime2FirstEvent,meanDwellTime,imh] = newRastergram(ah,traces2_allfr,ttb,fps,aligntype,delete,sorttype,sortState,color1)


color2=[0.9 0.9 0.9];
maxState = max(max(traces2_allfr));
% if  maxState == 1
%     color1=flip([0.28 	0.34 	0.65; 1 1 1]); % dark-blue intron-def   
% elseif maxState <=6
%     color1=[1 1 1; 0.28 	0.34 	0.65; 1 0 0; 0 0.4 0; 1 0 1; 0 1 0; 1 0.4 0.4];
% end
% set width of raster lines
width = 10;

% if filter == 1, filters out events > minfrm and events < maxfrm 
%in both data sets
filter=0; 
minfrm=10;
maxfrm=10000;

% figure formating (JEB)
set(0,'defaultlinelinewidth',0.3);

% if align == 1, aligns to first event in data set 1
if nargin < 4 || isempty(aligntype)
    aligntype = 0;
end
% if sortinv == 1, sorts data set
if nargin < 6 || isempty(sorttype)
    sorttype = 0;
end
% if delete == 1, deletes empty AOIs for data set 1(no events) 
if nargin < 5 || isempty(delete)
    delete = 0;
end

sortinv = (sorttype >0);
align = (aligntype >0);

if sorttype == 1
   aligntype = 1;     % aligns to first binding event
elseif sorttype == 2  
    aligntype = 2; % aligns to last departure event
elseif sorttype == 3
	align = 0;    % no alignment, this wouldn't make sense
end

% get initial # of AOIs
nAOIs = size(traces2_allfr,1);
% reduce number of AOIs that are displayed (JEB)
% nAOIs=100
% initialize evi vector
initial_nAOIs=nAOIs;
evi=1:nAOIs;

% initialize shift  == 0
shift = 0;
% cia = obj.createIntervals(traces2_allfr,ttb);

            
temp=traces2_allfr;
temp(traces2_allfr >= sortState) = 1;
temp(traces2_allfr < sortState) = 0;
cia = createIntervals(temp,ttb);
[tempCev,evi] = createCev(filter,cia,evi,sortinv,nAOIs,sorttype);
cev=tempCev.*sortState;
logik=[];
logik=(cia(:,1)==-2);
Time2FirstEvent=cia(logik,5); 
time_empty=max(Time2FirstEvent)*ones(length(Time2FirstEvent),1);
t=(Time2FirstEvent~=time_empty);
time21st=Time2FirstEvent(t);
meanTime2FirstEvent = mean(time21st);

logik=[];
logik=cia(:,1)==1;
DwellTime=cia(logik,5);
meanDwellTime = mean(DwellTime);

for i=1:max(max(traces2_allfr))
    if sortState~= i
        temp=traces2_allfr;
        temp(traces2_allfr == i) = 1;
        temp(traces2_allfr ~= i) = 0;
        cia2 = createIntervals(temp,ttb);
        [tempCev,evi] = createCev(filter,cia2,evi,sortinv,nAOIs,0);
        cev=cev+tempCev.*i;
    end
end


[cev,evi] = alignAndDeleteCev(delete,align,aligntype,cev,evi,sortState);
evi=flip(evi);
imh = plotRastergramFromCev(ah,initial_nAOIs,nAOIs,shift,fps,width,cev,color1,color2);
