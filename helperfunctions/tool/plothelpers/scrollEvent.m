classdef scrollEvent < handle
    properties
        src
        VerticalScrollCount
    end
    
    events
        scroll
    end
    methods
        function obj = scrollEvent(figure)
            set(figure, 'WindowScrollWheelFcn', @obj.scrollFun);
        end
        function scrollFun(obj,src,evt)
            obj.src=src;
            obj.VerticalScrollCount=evt.VerticalScrollCount;
            notify(obj,'scroll')
        end
    end
end
    
    