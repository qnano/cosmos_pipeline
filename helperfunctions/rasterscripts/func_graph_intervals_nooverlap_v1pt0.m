function graph=func_graph_intervals_nooverlap_v1pt0(h,file_path,header_path,min_frame,max_frame,spot_number,nodata,nopeak,peak,silent)
%file_path= file path
%the size of the grapgh will be [spot_number x (max_frame-min_frame)]
%if the area requested is bigger than the file, you will get a border of
%unused space.
%header_path= the path to a header file to use for the frames to time change. 
                %this should be a matlab .mat file
                %NOTE: IF you want the x axis to be frames, set header_path
                %to be empty. If you want time, set a header path.
%spot_number= number of spots to graph
%nodata= rgb matrix (1x3) of values 0-1. This will represent no data.
%nopeak= rgb matrix (1x3) of values 0-1. This will represent no peak.
%peak= rgb matrix (1x3) of values 0-1. This will represent a peak.



%example:
%fpath='C:\Users\Alex\Documents\MATLAB\hoskins100409_cy5tmp_600uW_intervals.dat'
%hpath='header.mat'
%out=func_graph_intervals(fpath,hpath,1223,2403,49,[0 0 0],[0 0 1], [0 1 0],0);
%this gives the entire file.

%silent=boolean value, if 1 then this function will not plot the figure.

%[name path]=uigetfile('D:\matlab\images','pick an interval file to load');
%if path==0 %do not proceed if uigetfile is cancelled
%    'cancelled loading file!'
%    return;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%V.1.1  Alex O   3/1/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now uses non consistent time information and graphs as rectangles instead
%of as a matrix
%V.2  Alex O   3/2/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now hard sets the correct size limit to reduce graphing errors
%V.2.1  Alex O   3/3/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now has inputs to set colors for events and choose between time or frames
%V.2.2  Alex O   3/8/10


loaded=load(file_path,'-mat');
intervals=loaded.Intervals.CumulativeIntervalArray;
%load header file and get frames to time info




graph=zeros(spot_number,max_frame-min_frame+1)-4;
spot=[];
spot_counter=0;
frame_counter=0;
for i=1:size(intervals,1)
   if intervals(i,1)== -2 || intervals(i,1)== -3 %start of new spot
       spot_counter=spot_counter+1;
       frame_counter=0;
       if spot_counter>spot_number %graphed the number of spots wanted
          break; 
       end
       %If not started at beginning will fill in gap here
       if intervals(i,2)>min(intervals(:,2)) && min_frame<intervals(i,2)
           gap=intervals(i,2)-min_frame; %the gap to fill (if min frame was the min of intervals(:,2)
           frame_diff=min_frame-min(intervals(:,2));%additional difference if min_frame is not = to intervals(:,2)
           if frame_diff>0
               frame_diff=0;
           end
           graph(spot_counter,1:gap+frame_diff)=-4;
           frame_counter=frame_counter+gap+frame_diff;
       end
   end
   %will fill in gaps between intervals here
   if i>1 %not first element
       gap=intervals(i,2)-intervals(i-1,3);
       if gap>1 %will now fill gap with -4
           graph(spot_counter,frame_counter:frame_counter+gap-1)=-4;
           frame_counter=frame_counter+gap-1;
       end
   end
   %now will set rgb values of pixels at (x,y)=(spot_counter,frame_counter)
   %set rgb values depending on spot
   for j=1:intervals(i,3)-intervals(i,2)+1 %for each frame of spot
       if j+intervals(i,2)-1>=min_frame %not less than min frame
           if j+intervals(i,2)-1<=max_frame %not greater than max frame
               frame_counter=frame_counter+1;
               graph(spot_counter,frame_counter)=intervals(i,1);
               %[spot_counter frame_counter  j intervals(i,1) j+intervals(i,2)-1]
           end
       end
   end
end
%now to replace all -2, 0, +2 with fours, 
graph(find(graph==-2))=4;
graph(find(graph==0))=4;
graph(find(graph==2))=4;
%now to replace all -3, +1, +3 with eights
graph(find(graph==-3))=16;
graph(find(graph==1))=16;
graph(find(graph==3))=16;
%Fixes graph to reflect incorrect requests
frame_diff=min(intervals(:,2))-min_frame;
if frame_diff>0 %requested min frame is lower than actual
    graph=[zeros(spot_number,frame_diff)-4 graph];
    graph(:,end-frame_diff+1:end)=[];
end
%in case of any errors, make it conform to the size limit
graph=graph(1:spot_number,1:max_frame-min_frame+1);
%if silent ==1 don;t graph anything, just return the matrix
if silent==1
   return; 
end
%get matrix of time info if header path exists
time2frames=[];
if ~isempty(header_path)
    header=load(header_path,'-mat');
    time2frames=header.vid.ttb(min_frame:max_frame)/1000;
end
%now to graph the results
func_graph_intervals_time_graph(h,graph,time2frames,-4,4,100,150,16,nodata,nopeak,[0 0 0],[0 0 0],peak);%100 and 150 shouldnt occur




