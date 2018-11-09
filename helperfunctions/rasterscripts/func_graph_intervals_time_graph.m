function func_graph_intervals_time_graph(h,graph,time2frames,sneg,s0,s1,s2,s3,cneg,c0,c1,c2,c3)
%This function plots the matrix graph as rectangles using non consistent
%time information from the header file of the experiment.
%graph=the matrix to be graphed [m by n]
%time2frames=time from first to last frame of graph, time2frames(i)=graph(1,i) [n by 1]
%
%File 1:    File2:  Result:(0=no peak, 1=peak -1=no data)
%-1         -1      -1 = red
%-1         0       -1
%-1         1       -1
%0          -1      -1
%0          0       0  = white    
%0          1       1  = blue
%1          -1      -1 
%1          0       2  = yellow
%1          1       3  = green
%sneg=-1, what value this is in the matrix graph
%s0=0,what value this is in the matrix graph
%s1=1,what value this is in the matrix graph
%s2=2,what value this is in the matrix graph
%s3=3,what value this is in the matrix graph
%cneg= rgb color matrix of -1
%c0= rgb color matrix of 0
%c1= rgb color matrix of 1
%c2= rgb color matrix of 2
%c3= rgb color matrix of 3

%For single files this will be only -1, 0 and 3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%V.1.0  Alex O   3/2/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now supports custom colors and frames or time graphs
%V.1.1  Alex O   3/8/10



%construct the graph
% figure;
for i=1:size(graph,1) %for each row of graph
    curr=graph(i,1);
    curr_length=1;
    for j=2:size(graph,2) %for each element of each row of graph
        if graph(i,j)==curr %current value still in place
            curr_length=curr_length+1;
        else %different value, graph current
            %graph color depending on value of curr
            if ~isempty(time2frames)
                xpos=time2frames(j-curr_length);
                width=time2frames(j)-time2frames(j-curr_length);
            else
                xpos=j-curr_length;
                width=curr_length;
            end
            if curr==sneg %-1 No data
                rectangle('Parent',h,'Position',[xpos,i-1,width,1],'EdgeColor',cneg,'FaceColor',cneg,'hittest','off');
            elseif curr==s0 %0 No peaks
                rectangle('Parent',h,'Position',[xpos,i-1,width,1],'EdgeColor',c0,'FaceColor',c0,'hittest','off');   
            elseif curr==s1 %1 %Peak on file2
                rectangle('Parent',h,'Position',[xpos,i-1,width,1],'EdgeColor',c1,'FaceColor',c1,'hittest','off');                
            elseif curr==s2% 2 %peak on file1
                rectangle('Parent',h,'Position',[xpos,i-1,width,1],'EdgeColor',c2,'FaceColor',c2,'hittest','off');
            elseif curr==s3% 3 %peak on both files (or just peak if only 1 file)
                rectangle('Parent',h,'Position',[xpos,i-1,width,1],'EdgeColor',c3,'FaceColor',c3,'hittest','off');
            else
                'this is an error of some kind!'
                curr
            end
            curr=graph(i,j);
            curr_length=1;
        end
    end
    %now to graph the last interval for each row
    %graph color depending on value of curr
    if ~isempty(time2frames)
        xpos=time2frames(j-curr_length);
        width=time2frames(j)-time2frames(j-curr_length);
    else
        xpos=j-curr_length;
        width=curr_length;
    end
    if curr==sneg %-1 No data
        rectangle('Parent',h,'Position',[xpos,i-1,width,1],'EdgeColor',cneg,'FaceColor',cneg,'hittest','off');
    elseif curr==s0 %0 No peaks
        rectangle('Parent',h,'Position',[xpos,i-1,width,1],'EdgeColor',c0,'FaceColor',c0,'hittest','off');   
    elseif curr==s1 %1 %Peak on file2
        rectangle('Parent',h,'Position',[xpos,i-1,width,1],'EdgeColor',c1,'FaceColor',c1,'hittest','off');                
    elseif curr==s2% 2 %peak on file1
        rectangle('Parent',h,'Position',[xpos,i-1,width,1],'EdgeColor',c2,'FaceColor',c2,'hittest','off');
    elseif curr==s3% 3 %peak on both files (or just peak if only 1 file)
        rectangle('Parent',h,'Position',[xpos,i-1,width,1],'EdgeColor',c3,'FaceColor',c3,'hittest','off');
    else
        'this is an error of some kind!'
        curr
    end
end
axis ij;
ylabel('Spots')
if ~isempty(time2frames)
    axis([time2frames(1) time2frames(end) 0 size(graph,1)]);
    xlabel('Time [s]')
else
    axis([1 size(graph,2) 0 size(graph,1)]);
    xlabel('Frames')
end

