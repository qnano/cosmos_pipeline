classdef analyzeMultiplexCompWrapperHMM < handle
   properties
        analyseDynamics
   end
   
    events
       ParamChange
   end
   
   methods
       
       function whoAmI(obj,src,~)
            basevars = evalin('base','whos');
            testClassvars = basevars(strcmp({basevars.class},class(obj)));

            for i = 1:length(testClassvars)
                if(eq(evalin('base',testClassvars(i).name),obj))
                    obj.name =testClassvars(i).name;
                end
            end
       end
       
        function obj = analyzeMultiplexCompWrapperHMM(h,gdlObj,gaObj,gdcObj,gsdObj,intergrationTime,colocalized)
            topTab = uitab(h, 'Title', 'Analyze');
            htabgroup = uitabgroup(topTab);
            menuhandle = uimenu('Label','Analyze');
            detectionFrame = 1;
            startFrame = 2;
            obj.analyseDynamics{1} = guiHMMAnalyse(htabgroup,gdlObj,gaObj,gdcObj,gsdObj,intergrationTime,menuhandle,'Component 1',startFrame,detectionFrame,colocalized,false);
            obj.analyseDynamics{2} = guiHMMAnalyse(htabgroup,gdlObj,gaObj,gdcObj,gsdObj,intergrationTime,menuhandle,'Component 2',startFrame,detectionFrame,colocalized,true);
     
            uimenu('Parent',menuhandle,'Label','Select Subregion','Callback',@obj.selectSubRegion);
            uimenu('Parent',menuhandle,'Label','Threshold Settings','Callback',@obj.setParamsGeneral);
            uimenu('Parent',menuhandle,'Label','Analyze','Callback',@obj.analyze);
            uimenu('Parent',menuhandle,'Label','Rastergram Settings','Callback',@obj.setRastergramSettings)
            uimenu('Parent',menuhandle,'Label','Calculate Rastergram','Callback',@obj.runRastergram);
            uimenu('Parent',menuhandle,'Label','Bootstrap Rastergram','Callback',@obj.bootStrapRastergramMenu);

        end   
        
        function setParamsGeneral(obj,~,~)
             obj.analyseDynamics{1}.setParamsGeneral;
             
             for i=2:size(obj.analyseDynamics,2)
                 obj.analyseDynamics{i}.THPmax=obj.analyseDynamics{1}.THPmax;
                 obj.analyseDynamics{i}.THP=obj.analyseDynamics{1}.THP;
                 obj.analyseDynamics{i}.THXmax=obj.analyseDynamics{1}.THXmax;
                 obj.analyseDynamics{i}.THX=obj.analyseDynamics{1}.THX;
                 obj.analyseDynamics{i}.THBGmax=obj.analyseDynamics{1}.THBGmax;
                 obj.analyseDynamics{i}.THBG=obj.analyseDynamics{1}.THBG;
                 obj.analyseDynamics{i}.THSmax=obj.analyseDynamics{1}.THSmax;
                 obj.analyseDynamics{i}.THS=obj.analyseDynamics{1}.THS;
                 obj.analyseDynamics{i}.thresholdChange;
             end
        end
        
       function setRastergramSettings(obj,~,~)
             obj.analyseDynamics{1}.setRastergramSettings;
             for i=2:size(obj.analyseDynamics,2)
                obj.analyseDynamics{i}.intergrationTime=obj.analyseDynamics{1}.intergrationTime;
                obj.analyseDynamics{i}.rastergramStartframe=obj.analyseDynamics{1}.rastergramStartframe;
                obj.analyseDynamics{i}.numberOfSpotsForRastergram=obj.analyseDynamics{1}.numberOfSpotsForRastergram;
             end
        end
        
        function bootStrapRastergramMenu(obj,~,~)
             obj.analyseDynamics{1}.bootStrapRastergramMenu;
             for i=2:size(obj.analyseDynamics,2)
                 obj.analyseDynamics{i}.bootStrapRastergram(obj.analyseDynamics{1}.Nbootstraps);
             end
        end
             
        function selectSubRegion(obj,C,~)
             for i=1:size(obj.analyseDynamics,2)
                 obj.analyseDynamics{i}.selectSubRegion(C);
             end
        end
                        
        
        function analyze(obj,src,~)
            for i=1:size(obj.analyseDynamics,2)
                obj.analyseDynamics{i}.analyze;
            end
        end
        
        function runRastergram(obj,src,~)
            for i=1:size(obj.analyseDynamics,2)
                obj.analyseDynamics{i}.runRastergram;
            end
        end
   end  
end