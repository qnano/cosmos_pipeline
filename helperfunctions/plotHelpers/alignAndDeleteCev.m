function [cev,evi] = alignAndDeleteCev(delete,align,aligntype,cev,evi,state)
    if nargin < 5 || isempty(state)
        state = -1;
    end
    %% delete empty events
    if delete ==1
        indx=[];
        for i=1:size(cev,1) 
%             if state == -1
                ind=find(cev(i,:)>0,1,'first'); % find rows with entries == 1
%             else
%                 ind=find(cev(i,:)==state,1,'first'); % find rows with entries == 1
%             end
            if ind~=0                        % if row has an entry == 1
                indx=[indx i];               % retain the row number
            end
        end
        cev=cev(indx,:);                         % keep only the rows with entries == 1
        evi=evi(indx);     
    end
      nAOIs=size(cev,1);                       % set number of AOIs to new reduced cev size
    %% align events
    if align == 1
        if aligntype == 1
            indx_f=[];
            for i=1:size(cev,1)
                if state == -1
                    indf=find(cev(i,:)>0,1,'first'); % find rows with column entries == 1
                else
                    indf=find(cev(i,:)==state,1,'first'); % find rows with column entries == 1
                end
                if indf~=0                        % if a row has an entry == 1 (an event)
                    indx_f=[indx_f indf];         % retain the first column number == 1    
                else
                    indx_f=[indx_f 1];            % else set column number == 1 
                end
            end
            max_f=max(indx_f);   % find max column number
            cev_new=[];          % again create new cumulative interval
            for i=1:size(cev,1)  % add new columns (== 10) such that all events align at first binding event
                % e.g. cev(1,:)= 0 0 0 0 1 1 1 0  0  0 
                %      index1=   1 2 3 4 5 7 8 9 10 11
                %      cev(2,:)= 0 0 0 1 1 1 1 1  0  0 
                %      index2=   1 2 3 4 5 7 8 9 10 11
                % max_f = 11, max_f-indx_f(1) = 11-5 = 6; max_f-indx_f(2) = 11-4 = 7; indx_f(1)= 5; indx_f(2)= 4 
                % Therefore: cev_new(1,:) = 10 10 10 10 10 10 0 0 0  0   1  1  1  0  0  0 10 10 10 10 10
                %                            1  2  3  4  5  6 7 8 9 10  11 12 13 14 15 16 17 18 19 20 21 
                %            cev_new(2,:) = 10 10 10 10 10 10 10 0 0  0  1  1  1  1  1  0  0  10 10 10 10
                %                            1  2  3  4  5  6  7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
                % both are aligned at index 11 
                cev_new=[cev_new; 10*ones(1,max_f-indx_f(i)) cev(i,:) 10*ones(1,indx_f(i))];
            end
            cev=cev_new; % rename cev_new to cev
            shift=max_f; % record shift of indices
        elseif  aligntype == 2
            indx_l=[];
            for i=1:size(cev,1)
                if state == -1
                    indl=find(cev(i,:)>0,1,'last');  % find rows with column entries == 1
                else
                    indl=find(cev(i,:)==state,1,'last');  % find rows with column entries == 1
                end
                if indl~=0                        % if row has an entry == 1
                    indx_l=[indx_l indl];         % retain the last column number == 1
                else
                    indx_l=[indx_l 1];            % else set column number == 1
                end
            end
            max_l=max(indx_l); % find max column number
            cev_new=[]; % again create new cumulative interval
            for i=1:size(cev,1)
                % add new columns (== 10) such that all events align at last departure event
                % e.g. cev(1,:)= 0 0 0 0 1 1 1 0  0  0 
                %      index1=   1 2 3 4 5 7 8 9 10 11
                %      cev(1,:)= 0 0 0 1 1 1 1 1  0  0 
                %      index2=   1 2 3 4 5 7 8 9 10 11
                % max_f = 11, max_f-indx_f(1) = 11-8 = 3; max_f-indx_f(2) = 11-9 = 2; indx_f(1)= 8; indx_f(2)= 9 
                % Therefore: cev_new(1,:) = 10 10 10  0  0  0  0 1 1  1  0  0  0 10 10 10 10 10 10 10 10 
                %                            1  2  3  4  5  6  7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 
                %            cev_new(2,:) = 10 10  0  0  0  1  1 1 1  1  0  0 10 10 10 10 10 10 10 10 10
                %                            1  2  3  4  5  6  7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
                % both are aligned at index 10
                cev_new=[cev_new; 10*ones(1,max_l-indx_l(i)) cev(i,:) 10*ones(1,indx_l(i))];
            end
            cev = cev_new; % rename cev_new to cev
            shift = max_l; % record shift of indices
         end
    end

    % invert_cev to sort rastergrams from bottom to top and not like
    % before frome top to bottom (JEB)
    cev = flipud(cev);
end