% This will make an Interval Data Structure from the traces (traces1_allfr and traces2_allfr) generated from 
% the output of camCal.m

vid=[0:1.5:898.5]; % this contains times of our frames (600 frames at 1.5s per frame)

traces1_allfr_bin=traces1_allfr;
traces1_allfr_bin(traces1_allfr_bin>0)=1;
traces1_allfr_bin(:,end)=[]; % because time refers to the start of the frame, need to do this so that each frame has start time and end time 
traces1_allfr_bins=num2str(traces1_allfr_bin); 
traces1_allfr_binsc=num2cell(traces1_allfr_bins,2); % cell array is needed for regexp to work
traces1_allfr_binsc=regexprep(traces1_allfr_binsc,' ',''); % this makes a searchable cell array where each cell is a binary trace in string format



Intervals.CumulativeIntervalArray=[]; %initialize array

for i=1:length(traces1_allfr_binsc)
% Key (Larry's flags):
% "1..." = -3
% "...1" = 3
% "0..." = -2
% "...0" = 2
% "1...0...1" = 0
% "0...1...0" = 1
exp=traces1_allfr_binsc{i,1};
[start_flag_3, end_flag_3]=regexp(exp,'^1+');
[start_flag_2, end_flag_2]=regexp(exp,'^0+');
[start_flag0, end_flag0]=regexp(exp,'(?<=1)0+(?=1)'); % only zeros flanked by ones
[start_flag1, end_flag1]=regexp(exp,'(?<=0)1+(?=0)'); % only ones flanked by zeros
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

save IDS1.mat Intervals
