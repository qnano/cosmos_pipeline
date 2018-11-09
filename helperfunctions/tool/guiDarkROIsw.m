classdef guiDarkROIsw < guiHMMAnalyse
    methods
        
        function obj = guiDarkROIsw(h,dataObject,alignmentObject,driftObject,spotDetectorObject,intergrationTime,menuhandle,menuName,startFrame,detectionFrame,colocalized,flip)
             if nargin < 8 || isempty(menuName)
                menuName = '';
            end
            
            if nargin <  6 || isempty(intergrationTime )
                intergrationTime = 1;
            end
            if nargin < 9 || isempty(startFrame)
                startFrame = 1;
            end
            if nargin < 10 || isempty(detectionFrame)
                detectionFrame = 1;
            end
            if nargin < 11 || isempty(colocalized)
                colocalized =1;
            end
            if nargin < 12 || isempty(flip)
                flip = 0;
            end
            if nargin < 7 || isempty(menuhandle)
                menuhandle=[];
            end
            if nargin < 8 || isempty(menuName)
            	menuName='Dark';
            else
                menuName= ['Dark' menuName];
            end
             if isa(spotDetectorObject,'spotDetectorMultiplexRNAWrapper')
                spotDetectorObject = spotDetectorObject.spotDetectorObjects{end};
             end
             
            obj@guiHMMAnalyse(h,dataObject,alignmentObject,driftObject,spotDetectorObject,intergrationTime,menuhandle,menuName,startFrame,detectionFrame,colocalized,flip);
            
            obj.classSingle = 'guiDarkROIsw';
            obj.classAdded = 'guiDarkROIswAdded';
        end
        
        function  analyze(obj,src,~)
            obj.clear;
            coords =   obj.spotDetectorObject.darkCoords;
            frames =   (obj.startFrame-1).*ones(size(coords,1),1);
            coords_1st=coords(frames==obj.startFrame-1,:);
     
             if ~obj.analyzeCoords1st(coords_1st)
                obj.plotImages;
            end
            obj.status=1;
        end
    end
end