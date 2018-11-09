classdef guiExperimentPoolingMultiplexRNAWrapper < analyzeMultiplexRNAWrapper
   properties

   end
   methods
        function obj = guiExperimentPoolingMultiplexRNAWrapper(varargin)
            colocalized = varargin{end};
            intergrationTime = varargin{end-1};
            MultiplexRNA = varargin{end-2};
            for i=1:MultiplexRNA
                sdmwObj.spotDetectorObjects{i}=[];
            end
            obj = obj@analyzeMultiplexRNAWrapper(varargin{1},[],[],[],sdmwObj,MultiplexRNA,intergrationTime,colocalized);
            obj.loadvars(varargin{2},false);
            for i=3:length(varargin)-3
                obj.addvars(varargin{i},false);
            end
        end  
   end
end