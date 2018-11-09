classdef dipTrackObj < handle
    %dipTrackObj  listen for updates to curslice value in userdata of dipimage figures
    
    properties
        slice; %old dipimage slice
    end
    
    methods
        
        function obj = dipTrackObj(slice)
            switch nargin
                case 0
                    obj.slice = 0;
                case 1
                    obj.slice = slice;
                otherwise
                    error('dipTrackObj:ConstructorImproperInputs',...
                        'dipTrackObj: 0 to 1 inputs required for constructor, dipTrackObj(figHandle,slice)')
            end
            %             addlistener(obj.figHandle,'UserData','PostSet',@obj.curslice);
        end
              
    
        %         function curslice(obj,h,evnt)
        %             %get userdata
        %             udata = get(evnt,'NewValue');
        %             %get current slice
        %             if isfield(udata,'curslice') && obj.slice ~= udata.curslice
        %                 obj.slice = udata.curslice; %update slice property
        %                 dipTrack(obj.figHandle,obj.slice); %update figure
        %             end
        %
        %         end
        
    end
    
end