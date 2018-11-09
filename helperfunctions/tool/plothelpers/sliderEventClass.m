classdef sliderEventClass < handle
    properties
        src
        VerticalScrollCount
        TH
        THmax
    end
    
    events
        sliderEvent
        sliderLimitsChange 
    end
    methods
        function obj = sliderEventClass(TH,THmax)
            obj.TH=TH;
            obj.THmax=THmax;
        end
        function changePosition(obj,TH)
            obj.TH=TH;
            notify(obj,'sliderEvent')
        end
        function changeLimits(obj,THmax)
            obj.THmax=THmax;
            notify(obj,'sliderLimitsChange')
        end
    end
end
    
    