classdef analyzeMultiplexRNAHMMWrapper < handle
   properties
        analyseDynamics
        objects
        addedvars=false;
        handels
   end
   
    events
       ParamChange
   end
   
   methods
       function deletetry(obj,handle)
            try
                delete(handle)
            catch
            end
       end
       
        function settry(obj,handle,var,prop)
            try
                set(handle,var,prop)
            catch
            end
        end
        
       function showTargetTrace(obj,bol)
           for i=1:size(obj.analyseDynamics,2)
                obj.analyseDynamics{i}.showTargetTrace(bol);
           end
       end
       
       function whoAmI(obj,~,~)
            basevars = evalin('base','whos');
            testClassvars = basevars(strcmp({basevars.class},class(obj)));
            
            for i = 1:length(testClassvars)
                if(eq(evalin('base',testClassvars(i).name),obj))
                    obj.name =testClassvars(i).name;
                end
            end
       end
        
       
        function obj = analyzeMultiplexRNAHMMWrapper(h,gdlObj,gaObj,gdcObj,sdmwObj,MultiplexRNA,intergrationTime,colocalized)
            if isa(h,'matlab.ui.container.Tab')
                htabgroup2 = uitabgroup(h); 
                htabgroup2.Visible='off';
            else
                htabgroup2 = h;
            end
            obj.objects.gdlObj = gdlObj;
            obj.objects.gaObj = gaObj;
            obj.objects.gdcObj = gdcObj;
            obj.objects.sdmwObj = sdmwObj;
                
            if isa(h,'matlab.ui.container.Tab')
                topTab=h;
            else
                topTab = uitab(h, 'Title', 'Analyze');
            end
            htabgroup = uitabgroup(topTab);
            obj.handels.parentMenu =  uimenu('Label','Analyze');
            
        
            for i=1:MultiplexRNA
                obj.analyseDynamics{i} = guiHMMAnalyse(htabgroup,obj.objects.gdlObj,obj.objects.gaObj,obj.objects.gdcObj,obj.objects.sdmwObj.spotDetectorObjects{i},intergrationTime,obj.handels.parentMenu,['RNA ' num2str(i)],MultiplexRNA+1,i,colocalized);
                obj.deletetry(obj.analyseDynamics{i}.handels.parentMenu)
            end
            
            obj.handels.childMenu(1) = uimenu('Parent',obj.handels.parentMenu,'Label','Select Subregion','Callback',@obj.selectSubRegion);
            obj.handels.childMenu(2) = uimenu('Parent',obj.handels.parentMenu,'Label','Threshold Settings','Callback',@obj.setParamsGeneral);
            obj.handels.childMenu(3) = uimenu('Parent',obj.handels.parentMenu,'Label','Rastergram Settings','Callback',@obj.setRastergramSettings)
            obj.handels.childCheckMenu(1) = uimenu('Parent',obj.handels.parentMenu,'Label','Analyze','Separator','on','Callback',@obj.analyze);
            obj.handels.childCheckMenu(2) = uimenu('Parent',obj.handels.parentMenu,'Label','Variational Bayes HMM Signal Fit','Callback',@obj.bayesSignalFitGMM);

            obj.handels.childMenu(10) = uimenu('Parent',obj.handels.parentMenu,'Label','Animation','Callback',@obj.exportAnimation);
            obj.handels.childMenu(5) = uimenu('Parent',obj.handels.parentMenu,'Label','Calculate Rastergram','Callback',@obj.runRastergram);
            obj.handels.childMenu(6) = uimenu('Parent',obj.handels.parentMenu,'Label','Bootstrap Rastergram','Callback',@obj.bootStrapRastergramMenu);
            obj.handels.childMenu(7) = uimenu('Parent',obj.handels.parentMenu,'Label','Load','Separator','on','Callback',@obj.loadvars);  
            obj.handels.childMenu(8) = uimenu('Parent',obj.handels.parentMenu,'Label','Add','Callback',@obj.addvars);              
            obj.handels.childMenu(9) = uimenu('Parent',obj.handels.parentMenu,'Label','Save','Callback',@obj.savevars); 

           for i=1:size(obj.analyseDynamics,2)
                obj.analyseDynamics{i}.handels.childCheckMenu(1)=obj.handels.childCheckMenu(1);
                obj.analyseDynamics{i}.handels.childCheckMenu(2)=obj.handels.childCheckMenu(2);
           end
        end  
        
        function exportAnimation(obj,~,~)
            Title = 'Export Animation Select RNA';
            Options.Resize = 'on';
            Options.Interpreter = 'tex';
            Options.CancelButton = 'on';
            Options.ApplyButton = 'off';
            Option.Dim = 1; % Horizontal dimension in fields

            Prompt = {};
            Formats = {};
            
            Prompt(1,:) = {'RNA number','index',[]};
            Formats(1,1).type = 'edit';
            Formats(1,1).format = 'integer';
            Formats(1,1).size = [-1 0];
            Formats(1,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.index = 1; 
            
            [settings,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);
            
            if ~Cancelled
                obj.analyseDynamics{settings.index}.exportAnimation;
            end
        
        end
         
        function bayesSignalFitGMM(obj,~,~)
            if ~obj.objects.sdmwObj.status
                error('Spots detector is not done.')
            else
                for i=1:size(obj.analyseDynamics,2)
                    obj.analyseDynamics{i}.bayesSignalFitGMM;
                end
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
        
        function bootStrapRastergram(obj,N)
            if nargin < 2 || isempty(N)
                N=10;
            end

             for i=1:size(obj.analyseDynamics,2)
                 obj.analyseDynamics{i}.bootStrapRastergram(N);
             end
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
        
        function selectSubRegion(obj,C,~)
            if nargin < 2
                C=[];
            end
            obj.analyseDynamics{1}.selectSubRegion(C);
            
            for i=1:size(obj.analyseDynamics,2)
                 obj.analyseDynamics{i}.subRegionSelection = obj.analyseDynamics{1}.subRegionSelection;
            end
        end
                        
        
        function analyze(obj,~,~)
            if ~obj.objects.sdmwObj.status
                error('Spots detector is not done.')
            else
                for i=1:size(obj.analyseDynamics,2)
                    obj.analyseDynamics{i}.analyze;
                end
            end
        end
        
        function runRastergram(obj,~,~)
            for i=1:size(obj.analyseDynamics,2)
                obj.analyseDynamics{i}.runRastergram;
            end
        end
        
        function vs = savevars(obj,src,~)
            

            vs.class = 'analyzeMultiplexRNAWrapper';
            for i=1:size(obj.analyseDynamics,2)
                vs.startFrame{i} = obj.analyseDynamics{i}.startFrame;
                vs.detectionFrame{i} = obj.analyseDynamics{i}.detectionFrame;
                vs.NFrames{i} = obj.analyseDynamics{i}.NFrames;
                vs.flipCameras{i} = obj.analyseDynamics{i}.flipCameras;
                vs.subRegionSelection{i} = obj.analyseDynamics{i}.subRegionSelection;
                vs.colocalized{i} = obj.analyseDynamics{i}.colocalized;
                vs.rastergramStartframe{i} = obj.analyseDynamics{i}.rastergramStartframe;
                vs.numberOfSpotsForRastergram{i} = obj.analyseDynamics{i}.numberOfSpotsForRastergram;
                vs.THP{i} = obj.analyseDynamics{i}.THP;
                vs.THX{i} = obj.analyseDynamics{i}.THX;
                vs.THBG{i} = obj.analyseDynamics{i}.THBG;
                vs.THS{i} = obj.analyseDynamics{i}.THS;
                vs.Nbootstraps{i} = obj.analyseDynamics{i}.Nbootstraps;
                vs.tagetTrace{i} = obj.analyseDynamics{i}.tagetTrace;
                vs.colorTarget{i} = obj.analyseDynamics{i}.colorTarget;
                vs.colorComplex{i}= obj.analyseDynamics{i}.colorComplex;


                %results
                vs.sortMap{i} = obj.analyseDynamics{i}.sortMap;
                vs.spotNrs{i} = obj.analyseDynamics{i}.spotNrs;
                vs.meanTime2FirstEvent{i} = obj.analyseDynamics{i}.meanTime2FirstEvent;
                vs.meanDwellTime{i} = obj.analyseDynamics{i}.meanDwellTime;
                vs.rawFitResultsCam1{i} = obj.analyseDynamics{i}.rawFitResultsCam1;
                vs.rawFitResultsCam2{i} = obj.analyseDynamics{i}.rawFitResultsCam2;
                vs.coordsCam1{i} = obj.analyseDynamics{i}.coordsCam1;
                vs.coordsCam2{i} = obj.analyseDynamics{i}.coordsCam2;
                vs.THPmax{i} = obj.analyseDynamics{i}.THPmax;
                vs.THXmax{i} = obj.analyseDynamics{i}.THXmax;
                vs.THSmax{i} = obj.analyseDynamics{i}.THSmax;
                vs.THBGmax{i} = obj.analyseDynamics{i}.THBGmax;


                vs.coord1{i} = obj.analyseDynamics{i}.coord1;
                vs.coord2{i} = obj.analyseDynamics{i}.coord2;
                vs.spotsIncluded{i} = obj.analyseDynamics{i}.spotsIncluded;
                vs.maskFilt2{i} = obj.analyseDynamics{i}.maskFilt2;
                vs.maskFilt1{i} = obj.analyseDynamics{i}.maskFilt1;
                vs.photons1I{i} = obj.analyseDynamics{i}.photons1I;
                vs.photons2I{i} = obj.analyseDynamics{i}.photons2I;
                vs.bg2I{i} = obj.analyseDynamics{i}.bg2I;
                vs.bg1I{i} = obj.analyseDynamics{i}.bg1I;
                vs.Sigma1I{i} = obj.analyseDynamics{i}.Sigma1I;
                vs.deltax{i} = obj.analyseDynamics{i}.deltax;
                vs.deltay{i} = obj.analyseDynamics{i}.deltay;
                vs.seq{i} = obj.analyseDynamics{i}.seq;
                vs.intergrationTime{i} = obj.analyseDynamics{i}.intergrationTime;
                vs.segmentedSpots{i} = obj.analyseDynamics{i}.segmentedSpots;
                vs.coordsComplex{i} = obj.analyseDynamics{i}.coordsComplex;
                vs.coordsTarget{i} = obj.analyseDynamics{i}.coordsTarget;
            end
            if nargin > 1 && isobject(src)
                prompt={'File name'};
                name = 'Save';
                defaultans = {['analyzeMultiplexRNAWrapper' date '.mat']};
                options.Interpreter = 'tex';
                answer = inputdlg(prompt,name,[1 40],defaultans,options);
                 save(answer{1},'vs','-v7.3')
%                 save(answer{1},'class','startFrame','detectionFrame','NFrames','flipCameras','subRegionSelection','colocalized','rastergramStartframe',...
%                     'numberOfSpotsForRastergram','THP','THX','THBG','THS','Nbootstraps','tagetTrace',...
%                     'sortMap','spotNrs','meanTime2FirstEvent','meanDwellTime','rawFitResultsCam1','rawFitResultsCam2','coordsCam1',...
%                     'coordsCam2','THPmax','THXmax','THSmax','THBGmax','coord1','coord2','spotsIncluded','maskFilt2','maskFilt1',...
%                     'photons1I','photons2I','bg2I','bg1I','Sigma1I','deltax','deltay','seq','intergrationTime','segmentedSpots',...
%                     'coordsComplex','coordsTarget')
            end
        end

        function loadvars(obj,src,updateplot)
            if ishandle(src)
                prompt={'File name'};
                name = 'Load';
                defaultans = {['analyzeMultiplexRNAWrapper' date '.mat']};
                options.Interpreter = 'tex';
                answer = inputdlg(prompt,name,[1 40],defaultans,options);
            elseif ischar(src)
                answer{1} = src;
            else
                answer{1} = [];  
            end
            if ~isempty(answer)
                if nargin < 3 || isempty(updateplot) || isobject(updateplot)
                    updateplot=true;
                end

                class = 'analyzeMultiplexRNAWrapper';
                if ischar(answer{1})
                    ws = load(answer{1});
                    ws=ws.vs;
                elseif isstruct(src)
                    ws=src;
                end
                if ~strcmpi(ws.class,class)
                    error('This is no analyzeMultiplexRNAWrapper class')
                else
                    for i=1:size(obj.analyseDynamics,2)
                        obj.analyseDynamics{i}.startFrame = ws.startFrame{i};
                        obj.analyseDynamics{i}.detectionFrame = ws.detectionFrame{i};
                        obj.analyseDynamics{i}.NFrames = ws.NFrames{i};
                        obj.analyseDynamics{i}.flipCameras = ws.flipCameras{i};
                        obj.analyseDynamics{i}.subRegionSelection = ws.subRegionSelection{i};
                        obj.analyseDynamics{i}.colocalized = ws.colocalized{i};
                        obj.analyseDynamics{i}.rastergramStartframe=ws.rastergramStartframe{i};
                        obj.analyseDynamics{i}.numberOfSpotsForRastergram=ws.numberOfSpotsForRastergram{i};
                        obj.analyseDynamics{i}.THP=ws.THP{i};
                        obj.analyseDynamics{i}.THX=ws.THX{i};
                        obj.analyseDynamics{i}.THBG=ws.THBG{i};
                        obj.analyseDynamics{i}.THS=ws.THS{i};
                        obj.analyseDynamics{i}.Nbootstraps=ws.Nbootstraps{i};
                        obj.analyseDynamics{i}.tagetTrace =ws.tagetTrace{i};
                        obj.analyseDynamics{i}.colorTarget =ws.colorTarget{i};
                        obj.analyseDynamics{i}.colorComplex = ws.colorComplex{i};

                        %results
                        obj.analyseDynamics{i}.sortMap=ws.sortMap{i};
                        obj.analyseDynamics{i}.spotNrs=ws.spotNrs{i};
                        obj.analyseDynamics{i}.meanTime2FirstEvent=ws.meanTime2FirstEvent{i};
                        obj.analyseDynamics{i}.meanDwellTime=ws.meanDwellTime{i};
                        obj.analyseDynamics{i}.rawFitResultsCam1=ws.rawFitResultsCam1{i};
                        obj.analyseDynamics{i}.rawFitResultsCam2=ws.rawFitResultsCam2{i};
                        obj.analyseDynamics{i}.coordsCam1=ws.coordsCam1{i};
                        obj.analyseDynamics{i}.coordsCam2=ws.coordsCam2{i};
                        obj.analyseDynamics{i}.THPmax=ws.THPmax{i};
                        obj.analyseDynamics{i}.THXmax=ws.THXmax{i};
                        obj.analyseDynamics{i}.THSmax=ws.THSmax{i};
                        obj.analyseDynamics{i}.THBGmax=ws.THBGmax{i};


                        obj.analyseDynamics{i}.coord1=ws.coord1{i};
                        obj.analyseDynamics{i}.coord2=ws.coord2{i};
                        obj.analyseDynamics{i}.spotsIncluded=ws.spotsIncluded{i};
                        obj.analyseDynamics{i}.maskFilt2=ws.maskFilt2{i};
                        obj.analyseDynamics{i}.maskFilt1=ws.maskFilt1{i};
                        obj.analyseDynamics{i}.photons1I=ws.photons1I{i};
                        obj.analyseDynamics{i}.photons2I=ws.photons2I{i};
                        obj.analyseDynamics{i}.bg2I=ws.bg2I{i};
                        obj.analyseDynamics{i}.bg1I=ws.bg1I{i};
                        obj.analyseDynamics{i}.Sigma1I=ws.Sigma1I{i};
                        obj.analyseDynamics{i}.deltax=ws.deltax{i};
                        obj.analyseDynamics{i}.deltay=ws.deltay{i};
                        obj.analyseDynamics{i}.seq=ws.seq{i};
                        obj.analyseDynamics{i}.intergrationTime=ws.intergrationTime{i};
                        obj.analyseDynamics{i}.segmentedSpots=ws.segmentedSpots{i};
                        obj.analyseDynamics{i}.coordsComplex=ws.coordsComplex{i};
                        obj.analyseDynamics{i}.coordsTarget=ws.coordsTarget{i};
                        obj.analyseDynamics{i}.handels.selectedSpotIndex = [];


                        maxNumberOfImages=size(obj.analyseDynamics{i}.spotsIncluded,1);
                        set(obj.analyseDynamics{i}.handels.spotSelectionSlider, 'SliderStep', [1/maxNumberOfImages 10/maxNumberOfImages]);       
                        set(obj.analyseDynamics{i}.handels.spotSelectionSlider, 'Max', maxNumberOfImages);       
                        set(obj.analyseDynamics{i}.handels.spotSelectionSlider, 'Value', round(maxNumberOfImages/2)); 

                        obj.analyseDynamics{i}.handels.Figure.UserData.sliderEventClass.changeLimits([obj.analyseDynamics{i}.THPmax; obj.analyseDynamics{i}.THXmax; obj.analyseDynamics{i}.THBGmax; obj.analyseDynamics{i}.THSmax]);
                        obj.analyseDynamics{i}.handels.Figure.UserData.sliderEventClass.changePosition([obj.analyseDynamics{i}.THP/(obj.analyseDynamics{i}.THPmax-1)*get(obj.analyseDynamics{i}.handels.intensityTHSlider,'Max');...
                        obj.analyseDynamics{i}.THX/(obj.analyseDynamics{i}.THXmax-1)*get(obj.analyseDynamics{i}.handels.positionTHSlider,'Max');...
                        obj.analyseDynamics{i}.THBG/(obj.analyseDynamics{i}.THBGmax-1)*get(obj.analyseDynamics{i}.handels.backgroundRatioTHSlider,'Max');...
                        obj.analyseDynamics{i}.THS/(obj.analyseDynamics{i}.THSmax-1)*get(obj.analyseDynamics{i}.handels.sigmaTHSlider,'Max')]);
                        obj.analyseDynamics{i}.draw3Dmovie=updateplot;
                            obj.analyseDynamics{i}.plotImages;
                    end
                end   
            end
        end
        
        function addvars(obj,src,updateplot)
            if ishandle(src)
                prompt={'File name'};
                name = 'Load';
                defaultans = {['analyzeMultiplexRNAWrapper' date '.mat']};
                options.Interpreter = 'tex';
                answer = inputdlg(prompt,name,[1 40],defaultans,options);
            elseif ischar(src)
                answer{1} = src;
            else
                answer{1} = [];  
            end
            if ~isempty(answer)
                if nargin < 3 || isempty(updateplot) || isobject(updateplot)
                    updateplot=true;
                end

                class = 'analyzeMultiplexRNAWrapper';
                if ischar(answer{1})
                    ws = load(answer{1});
                    ws=ws.vs;
                elseif isstruct(src)
                    ws=src;
                end
                if ~strcmpi(ws.class,class)
                    error('This is no analyzeMultiplexRNAWrapper class')
                else
                    for i=1:size(obj.analyseDynamics,2)
                        clear vs;
                        vs.class = 'guiHMMAnalyse';
                        vs.rawFitResultsCam1  = ws.rawFitResultsCam1{i};
                        vs.rawFitResultsCam2  = ws.rawFitResultsCam2{i};
                        vs.coordsCam1  = ws.coordsCam1{i};
                        vs.coordsCam2  = ws.coordsCam2{i};



                        vs.coord1  = ws.coord1{i};
                        vs.coord2  = ws.coord2{i};
                        vs.spotsIncluded  = ws.spotsIncluded{i};
                        vs.maskFilt2  = ws.maskFilt2{i};
                        vs.maskFilt1  = ws.maskFilt1{i};
                        vs.photons1I  = ws.photons1I{i};
                        vs.photons2I  = ws.photons2I{i};
                        vs.bg2I  = ws.bg2I{i};
                        vs.bg1I  = ws.bg1I{i};
                        vs.Sigma1I  = ws.Sigma1I{i};
                        vs.deltax  = ws.deltax{i};
                        vs.deltay  = ws.deltay{i};
                        vs.seq  = ws.seq{i};

                        vs.coordsComplex = ws.coordsComplex{i};
                        vs.coordsTarget = ws.coordsTarget{i};
                        obj.analyseDynamics{i}.addvars(vs,false);
                        
                        obj.settry(obj.handels.childMenu(5),'Separator','on')
                        obj.deletetry(obj.handels.childMenu(1))
                        obj.deletetry(obj.handels.childMenu(4))
                    

                    end
                end   
            end
        end

   end  
end