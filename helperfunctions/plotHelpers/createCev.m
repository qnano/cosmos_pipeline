function [cev,evi] = createCev(filter,cia,evi,sortinv,nAOIs,sorttype)
r=[];

%% set events that are too short or long to -1 to ignore
if filter==1
    % filter out events <= minfrm
    [r,~]=find(cia(:,1)==-3 & cia(:,4)<= minfrm);
    cia(r,1)=-1;  % sets indicator to -1
    [r,~]=find(cia(:,1)==1 & cia(:,4)<= minfrm);
    cia(r,1)=-1;  % sets indicator to -1
    [r,~]=find(cia(:,1)==3 & cia(:,4)<= minfrm);
    cia(r,1)=-1;  % sets indicator to -1
    
    % filter out events >= maxfrm
    [r,~]=find(cia(:,1)==-3 & cia(:,4)>= maxfrm);
    cia(r,1)=-1; % sets indicator to -1
    [r,~]=find(cia(:,1)==1 & cia(:,4)>= maxfrm);
    cia(r,1)=-1; % sets indicator to -1
    [r,~]=find(cia(:,1)==3 & cia(:,4)>= maxfrm);
    cia(r,1)=-1; % sets indicator to -1
end

%% sort events by duration
if sortinv == 1
     duration=zeros(1,nAOIs);
     firstbind=zeros(1,nAOIs);
     for i=1:nAOIs
        % find ith AOI events
        ia=cia(cia(:,7)==i,:);
        if sorttype == 1 || sorttype == 3 % find first event
            [r,~]=find(ia(:,1)==1 | ia(:,1)==-3 | ia(:,1)==3,1,'first');
        elseif sorttype == 2  % find last event
            [r,~]=find(ia(:,1)==1 | ia(:,1)==3,1,'last');
        end 
        if ~isempty(r)    % if there is a high event
            duration(i)= ia(r,4);     % get duration of that event
            firstbind(i) = ia(r,2);   % get first binding time of that event
        end
     end    
        
     if sorttype == 1 || sorttype == 2
         [ds,evi]=sort(duration);  % sort by duration
     elseif sorttype == 3     
         [ds,evi]=sort(firstbind); % sort by first binding time
         r=find(ds==0,1,'last');
         if ~isempty(r)
            evi = [evi(r+1:end) evi(1:r)]; % places all zero duration events (evi) at bottom of graph
         end
     end
     
%% create new cumulative intervals
    cia_new=[];
    for i=1:nAOIs
        % find ith AOI events
        ia=cia(cia(:,7)==i,:);
        cia_new=[cia_new; ia];  
    end   
    cia=cia_new;
end

%% create cumulative event matrix
% These cevs set all non-events (LJF: -2, 0 ,2; ADF: -1) equal to zero and all events (LJF -3,1,3) equal to 1
cev=[];
% loop over AOIs
for j=1:nAOIs
    % find jth AOI events
    ia=cia(cia(:,7)==evi(j),:);
    % loop over number of events
     for i=1:size(ia,1)
         % set all baselines == 0
         if (ia(i,1)==-2 || ia(i,1)==-1 || ia(i,1)==0 || ia(i,1)==2)
             cev(j,ia(i,2):ia(i,3)) = 0;
         %set all events == 1    
         elseif (ia(i,1)==-3 || ia(i,1)==1 || ia(i,1)==3)
              cev(j,ia(i,2):ia(i,3)) = 1;
         end
     end
end