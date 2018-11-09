classdef guiAnalyzeDynamics < handle
   properties
    handels
    PSFSigma=[];
    colorTarget;
    colorComplex;
    spotNrMap
    
    parent
    parentFigure
    
    hlo

    %settings
    startFrame
    detectionFrame
    NFrames
    flipCameras
    subRegionSelection
    colocalized
    rastergramStartframe
    numberOfSpotsForRastergram
    THP
    THX
    THBG
    THS
    Nbootstraps
    tagetTrace = true
    draw3Dmovie
    
    align
    delete
    sortType
    
    
    %results
    sortMap
    spotNrs
    meanTime2FirstEvent
    meanDwellTime
    rawFitResultsCam1
    rawFitResultsCam2
    coordsCam1
    coordsCam2
    THPmax
    THXmax
    THSmax
    THBGmax
    
    coord1
    coord2
    spotsIncluded
    maskFilt2
    maskFilt1
    photons1I
    photons2I
    bg2I
    bg1I
    Sigma1I
    deltax
    deltay
    seq
    intergrationTime
    segmentedSpots
    coordsComplex
    coordsTarget
    rasterGramOuput


    %object
    dataObject
    alignmentObject
    driftObject
    spotDetectorObject
    sliderEventClassObj
    
    status=0;
    addedvars=false;
    classSingle = 'guiAnalyzeDynamics'
    classAdded = 'guiAnalyzeDynamicsAdded'
    AddedTrackingInfo
   end
   
	events
       ParamChange
       timePointChanged
%        spotChanged
    end
    
    methods     
       function whoAmI(obj,~,~)
            basevars = evalin('base','whos');
            testClassvars = basevars(strcmp({basevars.class},class(obj)));
            
            for i = 1:length(testClassvars)
                if(eq(evalin('base',testClassvars(i).name),obj))
                    obj.name =testClassvars(i).name;
                end
            end
       end
       
       function selectSubRegion(obj,C,~)
            if size(C,1) == 2 && size(C,2) == 2
                obj.subRegionSelection = C;
            else
                h=dipshow(joinchannels('RGB',squeeze(obj.alignmentObject.rgbImage(:,:,1)),squeeze(obj.alignmentObject.rgbImage(:,:,2))),'lin');
                [~,obj.subRegionSelection] = dipcrop(h);
                close(h)
            end
       end
       
        function  setParamsGeneral(obj,~,~)
            prompt={'Enter max I','Enter current I','Enter max \Delta X','Enter current \Delta X','Enter max Photon x BG ','Enter current Photon x BG ','Enter max Sigma','Enter current Sigma'};
            name = 'Threshold Params';
               
            defaultans = {num2str(obj.THPmax),num2str(obj.THP), num2str(obj.THXmax),num2str(obj.THX), num2str(obj.THBGmax), num2str(obj.THBG),num2str(obj.THSmax), num2str(obj.THS)}; %1e3,10, 5
            options.Interpreter = 'tex';
            
            answer = inputdlg(prompt,name,[1 40],defaultans,options);
            if ~isempty(answer)
                obj.THPmax = max(str2num(answer{2}).*1.1,str2num(answer{1}));
                THP = str2num(answer{2});

                set(obj.handels.intensityTHSlider,'Value',THP/(obj.THPmax-1)*get(obj.handels.intensityTHSlider,'Max'));
                obj.THXmax = max(str2num(answer{4}).*1.1,str2num(answer{3}));
                THX = str2num(answer{4});
                set(obj.handels.positionTHSlider,'Value',THX/(obj.THXmax-1)*get(obj.handels.positionTHSlider,'Max'));
                obj.THBGmax = max(str2num(answer{6}).*1.1,str2num(answer{5}));
                THBG = str2num(answer{6});
                set(obj.handels.backgroundRatioTHSlider,'Value',THBG/(obj.THBGmax-1)*get(obj.handels.backgroundRatioTHSlider,'Max'));

                obj.THSmax = max(str2num(answer{8}).*1.1,str2num(answer{7}));
                THS = str2num(answer{8});
                set(obj.handels.sigmaTHSlider,'Value',THS/(obj.THSmax-1)*get(obj.handels.sigmaTHSlider,'Max'));

                obj.handels.Figure.UserData.sliderEventClass.changeLimits([obj.THPmax; obj.THXmax; obj.THBGmax; obj.THSmax]);
                obj.handels.Figure.UserData.sliderEventClass.changePosition([THP/(obj.THPmax-1)*get(obj.handels.intensityTHSlider,'Max');...
                    THX/(obj.THXmax-1)*get(obj.handels.positionTHSlider,'Max');...
                    THBG/(obj.THBGmax-1)*get(obj.handels.backgroundRatioTHSlider,'Max');...
                    THS/(obj.THSmax-1)*get(obj.handels.sigmaTHSlider,'Max')]);
            end
        end

       
        function obj = guiAnalyzeDynamics(h,dataObject,alignmentObject,driftObject,spotDetectorObject,intergrationTime,menuhandle,menuName,startFrame,detectionFrame,colocalized,flip)
            
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
            obj.parent = h;
            obj.status=0;
            obj.colocalized=colocalized;
            obj.flipCameras=flip;
            obj.intergrationTime=intergrationTime;
            obj.startFrame=startFrame;
            obj.detectionFrame=detectionFrame;
            
            obj.alignmentObject=alignmentObject;
            obj.driftObject=driftObject;
            obj.spotDetectorObject=spotDetectorObject; 
            
            obj.THPmax=1e3;
            obj.THXmax=10;
            obj.THSmax=10;
            obj.THBGmax=5;
            obj.rastergramStartframe=1;
            
            obj.align = false;
            obj.delete = false;
            obj.sortType = -1;
            
            obj.dataObject = dataObject;
            
            obj.handels.topTab = uitab(h, 'Title', ['Analyze ' menuName]);
            htabgroup = uitabgroup(obj.handels.topTab);
            obj.handels.traceEstimateTab = uitab(htabgroup, 'Title', 'Estimated Trace');
            obj.handels.selectedSpotTab = uitab(htabgroup, 'Title', 'Selected Spot');
            obj.handels.rastergramTap = uitab(htabgroup, 'Title', 'Rastergram');
            
            if nargin < 8 || isempty(menuName)
                menuName='Analyze';
            else
                menuName = ['Analyze ' menuName];
            end
            if nargin < 7 || isempty(menuhandle)
                obj.handels.parentMenu = uimenu('Label',menuName);
            else
                obj.handels.parentMenu = uimenu('Parent',menuhandle,'Label',menuName);
            end
            obj.handels.childMenu(1) = uimenu('Parent',obj.handels.parentMenu,'Label','Select Subregion','Callback',@obj.selectSubRegion);
            obj.handels.childMenu(2) = uimenu('Parent',obj.handels.parentMenu,'Label','Threshold Settings','Callback',@obj.setParamsGeneral);
            obj.handels.childMenu(3) = uimenu('Parent',obj.handels.parentMenu,'Label','Rastergram Settings','Callback',@obj.setRastergramSettings)
            obj.handels.childMenu(4) = uimenu('Parent',obj.handels.parentMenu,'Label','Analyze','Separator','on','Callback',@obj.analyze);
            obj.handels.childMenu(10) = uimenu('Parent',obj.handels.parentMenu,'Label','Animation','Callback',@obj.exportAnimation);
            obj.handels.childMenu(5) = uimenu('Parent',obj.handels.parentMenu,'Label','Calculate Rastergram','Callback',@obj.runRastergram);
            obj.handels.childMenu(6) = uimenu('Parent',obj.handels.parentMenu,'Label','Bootstrap Rastergram','Callback',@obj.bootStrapRastergramMenu);
            obj.handels.childMenu(11) = uimenu('Parent',obj.handels.parentMenu,'Label','Export Data','Callback',@obj.exportData);
            obj.handels.childMenu(7) = uimenu('Parent',obj.handels.parentMenu,'Label','Load','Separator','on','Callback',@obj.loadvars);  
            obj.handels.childMenu(8) = uimenu('Parent',obj.handels.parentMenu,'Label','Add','Callback',@obj.addvars);              
            obj.handels.childMenu(9) = uimenu('Parent',obj.handels.parentMenu,'Label','Save','Callback',@obj.savevars);
             
            obj.handels.timePointListener = addlistener(obj,'timePointChanged',@obj.plotTimePoint)
            
            slmin = 1;
            slmax = 100;
            obj.handels.spotSelectionSlider = uicontrol('Parent', obj.handels.traceEstimateTab,'Style','slider','Min',slmin,'Max',slmax,...
                            'SliderStep',[1 1]./(slmax-slmin),'Value',51,...
                            'Units','normalized','Position',[ 0.1 0.01 0.8 0.05]);
     
            slmin = 0;
            slmax = 100;
            obj.handels.intensityTHSlider = uicontrol('Parent',  obj.handels.traceEstimateTab,'Style','Slider','SliderStep',[1 1]./(slmax-slmin)/2,...
            'Units','normalized','Position',[0 0 0.01 0.8],'Min',slmin,'Max',slmax,...
                  'Value',51,'String','a');                           
            obj.handels.positionTHSlider = uicontrol('Parent',  obj.handels.traceEstimateTab,'Style','Slider','SliderStep',[1 1]./(slmax-slmin)/2,...
            'Units','normalized','Position',[0.01 0 0.01 0.8],'Min',slmin,'Max',slmax,...
                  'Value',51); 
            obj.handels.backgroundRatioTHSlider = uicontrol('Parent',  obj.handels.traceEstimateTab,'Style','Slider','SliderStep',[1 1]./(slmax-slmin)/2,...
            'Units','normalized','Position',[0.02 0 0.01 0.8],'Min',slmin,'Max',slmax,...
                  'Value',0); 
            obj.handels.sigmaTHSlider = uicontrol('Parent',  obj.handels.traceEstimateTab,'Style','Slider','SliderStep',[1 1]./(slmax-slmin)/2,...
            'Units','normalized','Position',[0.03 0 0.01 0.8],'Min',slmin,'Max',slmax,...
                  'Value',51); 
            
           
            obj.handels.settingsTextTH = uicontrol('Parent',  obj.handels.traceEstimateTab,'Style','text','Units','normalized','position',[ 0 0.8 0.05 0.1],'String','');
            obj.handels.spotIncludedCheckBox = uicontrol('Parent',  obj.handels.traceEstimateTab,'Style','checkbox',... 
                'String','Include trace','Units','normalized',... 
                'Value',0,'Position',[0 .9 0.1 0.1]);
            
            set(obj.handels.spotIncludedCheckBox,'Callback',@obj.includeSpot);       
            set(obj.handels.spotSelectionSlider,'Callback',@obj.changeSpot);
            
       
            obj.handels.Figure = getParentFigure(h);
            if ~isfield(obj.handels.Figure.UserData,'sliderEventClass')
                obj.handels.Figure.UserData.sliderEventClass = sliderEventClass(...
                    [get(obj.handels.intensityTHSlider,'Value'); get(obj.handels.positionTHSlider,'Value'); get(obj.handels.backgroundRatioTHSlider,'Value'); get(obj.handels.sigmaTHSlider,'Value')],...
                    [ obj.THPmax;obj.THXmax;obj.THBGmax; obj.THSmax])
            end
            obj.hlo = addlistener(obj.handels.Figure.UserData.sliderEventClass,'sliderEvent',@obj.thresholdChange)
            
            set(obj.handels.intensityTHSlider,'Callback',@obj.changePosition);
            set(obj.handels.positionTHSlider,'Callback',@obj.changePosition);
            set(obj.handels.backgroundRatioTHSlider,'Callback',@obj.changePosition);
            set(obj.handels.sigmaTHSlider,'Callback',@obj.changePosition);
             
        end
        
        function exportAnimation(obj,src,~)
            if obj.status == 1
                   
            %%%% SETTING DIALOG OPTIONS
            % Options.WindowStyle = 'modal';
            Title = 'Export Animation';
            Options.Resize = 'on';
            Options.Interpreter = 'tex';
            Options.CancelButton = 'on';
            Options.ApplyButton = 'off';
            Option.Dim = 1; % Horizontal dimension in fields

            Prompt = {};
            Formats = {};

            Prompt(1,:) = {'Spot IDs','spotIdString',[]};
            Formats(1,1).type = 'edit';
            Formats(1,1).format = 'text';
            Formats(1,1).size = [-1 0];
            Formats(1,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.spotIdString = '1:10'; % [pwd '/Data.tif'];
            
            Prompt(2,:) = {'\Delta T','dt',[]};
            Formats(2,1).type = 'edit';
            Formats(2,1).format = 'integer';
            Formats(2,1).size = [-1 0];
            Formats(2,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.dt = 1; %[pwd '/Dark.tif'];

            Prompt(3,:) = {'Frame rate','frameRate',[]};
            Formats(3,1).type = 'edit';
            Formats(3,1).format = 'integer';
            Formats(3,1).size = [-1 0];
            Formats(3,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.frameRate = 15; %[pwd '/Dark.tif'];

            Prompt(4,:) = {'Filename','exportAnimationPath',[]};
            Formats(4,1).type = 'edit';
            Formats(4,1).format = 'file';
            Formats(4,1).limits = [1 0];
            Formats(4,1).size = [-1 0];
            Formats(4,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.exportAnimationPath = [pwd '\exportAnimation' date '.avi']; % [pwd '/Data.tif'];           
            
            [settings,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);
         
                if ~Cancelled
                    Sample = settings.dt;
                    NFrame = obj.rawFitResultsCam1.Frame(end)+1;
                    frames2 = obj.rawFitResultsCam2.Frame(obj.maskFilt2,:);
                    frames1 = obj.rawFitResultsCam1.Frame(obj.maskFilt1,:);

                    visSpots = eval(settings.spotIdString);

                    h = figure('position', [100, 100, 800, 700])

                    diptruesize(h,200)
                    
                    % Prepare VideoWriter object
                    TimeStep=settings.frameRate;
                    movieFileName2 = settings.exportAnimationPath;
                    vidObj = VideoWriter(movieFileName2);
                    vidObj.FrameRate = round(TimeStep);

                    %open avi file
                    open(vidObj);
                    traces2_allfr = obj.getBinary(ones(size(obj.photons1I,1),1));

                    for jj=1:length(visSpots)
                    Nspot = visSpots(jj);
                        for ii = 0:Sample:NFrame-1

                            idx = find(frames2==ii);
                            A(:,:,ii+1) = obj.rawFitResultsCam2.ROIStack(:,:,idx(Nspot));
                            idx1 = find(frames1==ii);
                            A1(:,:,ii+1) = obj.rawFitResultsCam1.ROIStack(:,:,idx1(Nspot));
                            xxyy(ii+1,:) = obj.rawFitResultsCam2.Coord(idx(Nspot),:)-obj.rawFitResultsCam2.RoiStart(idx(Nspot),:);
                            xxyytarget(ii+1,:) = obj.coordsCam2(idx(Nspot),1:2)-obj.rawFitResultsCam2.RoiStart(idx(Nspot),:);

                            h1 = subplot(2,1,1)
                            plot(obj.photons1I(Nspot,:,1+obj.flipCameras),['-'],'color',obj.colorTarget)
                            hold on
                            plot(obj.photons1I(Nspot,:,2-obj.flipCameras),['-'],'color', obj.colorComplex)
                            plot(obj.bg1I(Nspot,:,2-obj.flipCameras)*(obj.PSFSigma*2)^2,'Color',[.7 .5 0])

                            plot(traces2_allfr(Nspot,:)*max(max(obj.photons1I(Nspot,:,2-obj.flipCameras)),max(obj.photons1I(Nspot,:,1+obj.flipCameras))),'-b')
                            ylabel('I [# Photons]')
                            xlabel(['Frame x ' num2str(obj.intergrationTime) ' [s]'])
                            axis tight
                            legend('Complex','Target','Background','Binary','Location','southoutside','Orientation','horizontal')

                            xlim([1, size(obj.photons1I(Nspot,:,2-obj.flipCameras),2)]);

                          if ~isempty(obj.seq)
                                b = obj.seq(Nspot,:);
                                %title(sprintf('spot ID = %d, Frame = %d, \\Delta X = %0.2g, I = %0.2g, bg=%0.2g',obj.handels.selectedSpotIndex,ii,b(ii),obj.photons1I(obj.handels.selectedSpotIndex,ii),obj.bg1I(obj.handels.selectedSpotIndex,ii,2)*(obj.PSFSigma*2)^2),'Parent',obj.handels.traceEstimateTabPlot);
                                h2= title(sprintf('spot ID = %d\n Frame = %d, \\Delta X = %0.2g, I = %0.2g, bg=%0.2g, S=%0.2g',...
                                    Nspot,ii+1,sqrt(sum((xxyy(ii+1,:)-xxyytarget(ii+1,:)).^2)),obj.photons1I(Nspot,ii+1,2-obj.flipCameras),...
                                    obj.bg1I(Nspot,ii+1,2-obj.flipCameras)*(obj.PSFSigma*2)^2,...
                                    obj.Sigma1I(Nspot,ii+1,2-obj.flipCameras)),...
                                    'Parent',h1);
                            end              

                            vline(ii,'-k','t')
                            hold off
                            hs = subplot(2,1,2);

                            axes1Position = get(gca, 'Position');
                            delete(hs);

                            % Position the logo in the upper right.
                            x1=0;
                            y1=0.01;
                            hAxis2 = axes('Position', [x1 y1 0.5 0.5]);
                            axis off; % Turn off tick marks, etc.
                            gimr = uint8(256.*cat(3,obj.colorComplex(1).*mat2gray(A(:,:,ii+1)'), obj.colorComplex(2).*mat2gray(A(:,:,ii+1)'), obj.colorComplex(3).*mat2gray(A(:,:,ii+1)'))); %flip(flip(A(:,:,ii+1),1),2));
                            
                            imshow(gimr);
                            hold on
                            plot(xxyy(ii+1,1)+1,xxyy(ii+1,2)+1,'h','MarkerFaceColor',obj.colorComplex,'MarkerEdgeColor','k','MarkerSize',8)
                            plot(xxyytarget(ii+1,1)+1,xxyytarget(ii+1,2)+1,'h','MarkerFaceColor',obj.colorTarget,'MarkerEdgeColor','k','MarkerSize',8)
                            title('Complex')
                            hAxis2 = axes('Position', [x1+0.5 y1 0.5 0.5]);

                            gimg = mat2gray(flip(flip(A1(:,:,ii+1),1),2)); 
                            gimg = uint8(gimg * 256);
                            gimg = cat(3,obj.colorTarget(1).*gimg,obj.colorTarget(2).*gimg,obj.colorTarget(3).*gimg);

                            rgbImage = gimr+gimg;


                            imshow(rgbImage);
                            hold on
                            plot(-10,-10,'h','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8)
                            plot(xxyy(ii+1,1)+1,xxyy(ii+1,2)+1,'h','MarkerFaceColor',obj.colorComplex,'MarkerEdgeColor','k','MarkerSize',8)
                            plot(xxyytarget(ii+1,1)+1,xxyytarget(ii+1,2)+1,'h','MarkerFaceColor',obj.colorTarget,'MarkerEdgeColor','k','MarkerSize',8)

                            legend('Position')
                            title('Target & Complex')
                            currFrame = getframe(h);
                            data = uint8(extend(currFrame.cdata,[194 194 3]));
                            currFrame.cdata = data;
                            writeVideo(vidObj,currFrame);
                            obj.deletetry(h2);
                        end
                    end
                     %close avi file
                close(vidObj);
                close(h);
                end
            end
        end
        
        function thresholdVec = getThreshold(obj,src,~)
            thresholdVec(1) = obj.THP;
            thresholdVec(2) = obj.THX;
            thresholdVec(3) = obj.THBG;
            thresholdVec(4) = obj.THS;
        end
        
        function setThreshold(obj,thresholdVec)           
                obj.THP = thresholdVec(1);
                obj.THX = thresholdVec(2);
                obj.THBG = thresholdVec(3);
                obj.THS = thresholdVec(4);
                
                obj.THPmax = max(obj.THP.*1.1,obj.THPmax);
                obj.THXmax = max(obj.THX.*1.1,obj.THXmax);
                obj.THBGmax = max(obj.THBG.*1.1,obj.THBGmax);
                obj.THSmax = max(obj.THS.*1.1,obj.THSmax);
                
                set(obj.handels.intensityTHSlider,'Value',obj.THP/(obj.THPmax-1)*get(obj.handels.intensityTHSlider,'Max'));
                set(obj.handels.positionTHSlider,'Value',obj.THX/(obj.THXmax-1)*get(obj.handels.positionTHSlider,'Max'));
                set(obj.handels.backgroundRatioTHSlider,'Value',obj.THBG/(obj.THBGmax-1)*get(obj.handels.backgroundRatioTHSlider,'Max'));
                set(obj.handels.sigmaTHSlider,'Value',obj.THS/(obj.THSmax-1)*get(obj.handels.sigmaTHSlider,'Max'));

                obj.handels.Figure.UserData.sliderEventClass.changeLimits([obj.THPmax; obj.THXmax; obj.THBGmax; obj.THSmax]);
                obj.handels.Figure.UserData.sliderEventClass.changePosition([...
                    obj.THP/(obj.THPmax-1)*get(obj.handels.intensityTHSlider,'Max');...
                    obj.THX/(obj.THXmax-1)*get(obj.handels.positionTHSlider,'Max');...
                    obj.THBG/(obj.THBGmax-1)*get(obj.handels.backgroundRatioTHSlider,'Max');...
                    obj.THS/(obj.THSmax-1)*get(obj.handels.sigmaTHSlider,'Max')]);
        end
        
        function setLimits(obj,src,ev)
            obj.THPmax = obj.handels.Figure.UserData.sliderEventClass.THmax(1);
            obj.THXmax = obj.handels.Figure.UserData.sliderEventClass.THmax(2);
            obj.THBGmax = obj.handels.Figure.UserData.sliderEventClass.THmax(3);
            obj.THSmax = obj.handels.Figure.UserData.sliderEventClass.THmax(4);
        end
        
        function changePosition(obj,src,ev) 
            TH(1) = get(obj.handels.intensityTHSlider,'Value')
            TH(2) = get(obj.handels.positionTHSlider,'Value')
            TH(3) = get(obj.handels.backgroundRatioTHSlider,'Value')
            TH(4) = get(obj.handels.sigmaTHSlider,'Value')
            
            obj.handels.Figure.UserData.sliderEventClass.changePosition(TH);
            obj.thresholdChange;
        end
        
        function includeSpot(obj,src,~)
           if ~isempty(obj.spotsIncluded) && ~isempty(obj.handels.selectedSpotIndex)
                obj.spotsIncluded(obj.handels.selectedSpotIndex) = get(obj.handels.spotIncludedCheckBox,'Value');
            end
            obj.updateBoxes;
        end
        
        function updateBoxes(obj)
            
             boxsize=10;
            
              if isfield(obj.handels,'analyzedTargets') % && ~isempty(obj.handels.analyzedTargets) && 2 <= size(obj.handels.analyzedTargets,2) && ishandle(obj.handels.analyzedTargets(1))
                    for i=1:length(obj.handels.analyzedTargets)
                        obj.deletetry(obj.handels.analyzedTargets(i))
                    end
              end
             
             if ~isempty(obj.spotsIncluded) && ~isempty(obj.handels.selectedSpotIndex)
                 if ~isempty(obj.coord1)  && isfield(obj.handels,'imageStackCam1')
                     obj.handels.analyzedTargets(1) = line(obj.coord1(obj.spotsIncluded,1)+1,obj.coord1(obj.spotsIncluded,2)+1,...
                    'Color','g','marker','o','LineStyle','none','Parent',obj.handels.imageStackCam1,'markersize',boxsize-1,'tag','pointsSM','buttondownfcn',@(src,eventdata)obj.selectMarker(src,eventdata,1));
                 end
                 if ~isempty(obj.coord2) && isfield(obj.handels,'imageStackCam2')
                    obj.handels.analyzedTargets(3) = line(obj.coord2(obj.spotsIncluded,1)+1,obj.coord2(obj.spotsIncluded,2)+1,...
                    'Color','g','marker','o','LineStyle','none','Parent',obj.handels.imageStackCam2,'markersize',boxsize-1,'tag','pointsSM','buttondownfcn',@(src,eventdata)obj.selectMarker(src,eventdata,2));                
                end

                if sum(~obj.spotsIncluded) > 0 
                    obj.handels.analyzedTargets(2) = line(obj.coord1(~obj.spotsIncluded,1)+1,obj.coord1(~obj.spotsIncluded,2)+1,...
                    'Color','r','marker','o','LineStyle','none','Parent',obj.handels.imageStackCam1,'markersize',boxsize-1,'tag','pointsSM','buttondownfcn',@(src,eventdata)obj.selectMarker(src,eventdata,1));
                    obj.handels.analyzedTargets(4) = line(obj.coord2(~obj.spotsIncluded,1)+1,obj.coord2(~obj.spotsIncluded,2)+1,...
                    'Color','r','marker','o','LineStyle','none','Parent',obj.handels.imageStackCam2,'markersize',boxsize-1,'tag','pointsSM','buttondownfcn',@(src,eventdata)obj.selectMarker(src,eventdata,2));
               end
             end

        end
        
        
        function changeSpot(obj,src,~)
            if  ~isempty(squeeze(obj.coord1))
                pos = round(get(src,'Value'));
                obj.handels.selectedSpotIndex = round((size(obj.coord1,1))/get(src,'Max')*pos);
%                 notify(obj,'spotChanged');
                obj.plotImages;
            end
        end
        
        function dataChange(obj,src,~)
            notify(obj,'ParamChange');
        end
      
        function  analyze(obj,src,~)
            if ~obj.dataObject.status
                error('No data loaded!')
            end
            if ~obj.alignmentObject.status 
                error('No alignment completed!')
            end
            if ~obj.driftObject.status 
                 error('No drift completed!')
            end
            if ~obj.spotDetectorObject.status
                error('Spots detector is not done.')
            end
            
            if isfield(obj.dataObject.fileStructure,'Data')
                obj.AddedTrackingInfo.files=[];
                obj.AddedTrackingInfo.files{1}  = obj.dataObject.fileStructure.Data;
            end
            
            obj.clear;
            coords =  obj.spotDetectorObject.rawInitialFitResultsCam1.Coord(obj.spotDetectorObject.maskPreFiltCam1 & obj.spotDetectorObject.maskFilt1,:);
            frames =  obj.spotDetectorObject.rawInitialFitResultsCam1.Frame(obj.spotDetectorObject.maskPreFiltCam1 & obj.spotDetectorObject.maskFilt1,:);
            coords_1st=coords(frames==obj.detectionFrame-1,:);

            if ~obj.analyzeCoords1st(coords_1st)
                obj.plotImages;
            end
            obj.status=1;
        end
        
        function changed = analyzeCoords1st(obj,coords_1st)
        
            tSegmentedSpots = transformPointsInverse(obj.alignmentObject.tformTotal,coords_1st);
            if ~isempty(obj.subRegionSelection)
                idx = obj.subRegionSelection(1,1) < tSegmentedSpots(:,1) &...
                obj.subRegionSelection(1,1)+obj.subRegionSelection(2,1) > tSegmentedSpots(:,1)&...
                obj.subRegionSelection(1,2) < tSegmentedSpots(:,2) &...
                obj.subRegionSelection(1,2)+obj.subRegionSelection(2,2) > tSegmentedSpots(:,2);
            else
                idx = true(size(tSegmentedSpots,1),1);
            end
            if sum(idx)>0
                maxNumberOfImages=sum(idx);
                obj.spotsIncluded = true(maxNumberOfImages,1);
                set(obj.handels.spotSelectionSlider, 'SliderStep', [max(0,min(1,1/maxNumberOfImages)) max(0,min(1,10/maxNumberOfImages))]);       
                set(obj.handels.spotSelectionSlider, 'Max', maxNumberOfImages);       
                set(obj.handels.spotSelectionSlider, 'Value', round(maxNumberOfImages/2)); 
                obj.segmentedSpots = coords_1st(idx,:);

                %%
                obj.coordsCam1=[];
                obj.coordsCam2=[];
                obj.colorTarget = obj.dataObject.getColor(1);
                obj.colorComplex = obj.dataObject.getColor(2);
                obj.intergrationTime = 1./obj.dataObject.getFramesPerSecond;

                for w=obj.detectionFrame:size(obj.dataObject.getData(1),3)
                    coords_nth = obj.segmentedSpots-repmat(obj.driftObject.drift(w,:)-obj.driftObject.drift(obj.detectionFrame,:),size(obj.segmentedSpots,1),1);
                    if w>=obj.startFrame
                        obj.coordsCam1=cat(1,obj.coordsCam1,[coords_nth (w-1).*ones(size(obj.segmentedSpots,1),1)]);
                        obj.coordsCam2=cat(1,obj.coordsCam2,[transformPointsInverse(obj.alignmentObject.tformTotal,coords_nth) (w-1).*ones(size(obj.segmentedSpots,1),1)]); 
                    end
                end

                % MLE Fit Intensities Cam1
                paramsFit = obj.spotDetectorObject.paramsFit;
                paramsFit.FitSigma=true;
                 coords1 = [obj.coordsCam1(:,1) obj.coordsCam1(:,2) obj.coordsCam1(:,3)];
                [ obj.rawFitResultsCam1 ] = fitBoxCenters( single(squeeze(obj.dataObject.getData(1))),coords1,paramsFit);

                 coords2 = [obj.coordsCam2(:,1) obj.coordsCam2(:,2) obj.coordsCam2(:,3)];
                [ obj.rawFitResultsCam2 ] = fitBoxCenters( single(squeeze(obj.dataObject.getData(2))),coords2,paramsFit);
                obj.PSFSigma=obj.spotDetectorObject.paramsFit.PSFSigma;

                obj.NFrames  =obj.rawFitResultsCam2.Frame(end)+2-obj.startFrame;
                obj.maskFilt2 = true(size(obj.rawFitResultsCam2.Photons,1),1);
                obj.maskFilt1 = true(size(obj.rawFitResultsCam1.Photons,1),1);

                photons1 = obj.rawFitResultsCam1.Photons(obj.maskFilt1,:);
                bg1 = obj.rawFitResultsCam1.Bg(obj.maskFilt1,:);
                bg2 = obj.rawFitResultsCam2.Bg(obj.maskFilt2,:);
                photons2 = obj.rawFitResultsCam2.Photons(obj.maskFilt2,:);

                obj.photons1I(:,:,1+obj.flipCameras) = reshape(photons1,[size(photons1,1)/obj.NFrames obj.NFrames]);
                obj.photons1I(:,:,2-obj.flipCameras) = reshape(photons2,[size(photons1,1)/obj.NFrames obj.NFrames]);

                sigma1 = min(obj.rawFitResultsCam1.Sigma(obj.maskFilt1,1),obj.THSmax-9*eps);
                sigma2 = min(obj.rawFitResultsCam2.Sigma(obj.maskFilt2,1),obj.THSmax-9*eps);
                obj.Sigma1I(:,:,1+obj.flipCameras) = reshape(sigma1,[size(sigma1,1)/obj.NFrames obj.NFrames]);
                obj.Sigma1I(:,:,2-obj.flipCameras) = reshape(sigma2,[size(sigma2,1)/obj.NFrames obj.NFrames]);

                obj.bg1I(:,:,2-obj.flipCameras) = reshape(bg2,[size(photons2,1)/obj.NFrames obj.NFrames]);
                obj.bg1I(:,:,1+obj.flipCameras) = reshape(bg1,[size(photons2,1)/obj.NFrames obj.NFrames]);

                nn=size(obj.rawFitResultsCam2.RoiStart(obj.maskFilt2,1),1);

                obj.coordsComplex = obj.rawFitResultsCam2.Coord(obj.maskFilt2,:)-obj.rawFitResultsCam2.RoiStart(obj.maskFilt2,:);
                obj.coordsTarget = obj.coordsCam2(obj.maskFilt2,1:2)-obj.rawFitResultsCam2.RoiStart(obj.maskFilt2,:);

                delta = obj.coordsComplex-obj.coordsTarget;
                obj.deltax = reshape(delta(obj.maskFilt2,1),[size(delta(obj.maskFilt2,1),1)/obj.NFrames obj.NFrames]);
                obj.deltay = reshape(delta(obj.maskFilt2,2),[size(delta(obj.maskFilt2,2),1)/obj.NFrames obj.NFrames]);
                obj.coordsComplex = reshape(obj.coordsComplex,[nn/obj.NFrames obj.NFrames 2]);
                obj.coordsTarget = reshape(obj.coordsTarget,[nn/obj.NFrames obj.NFrames 2]);

                obj.seq = sqrt(obj.deltax.^2+obj.deltay.^2);

                changed = obj.thresholdChange;
            else
                changed = true;
            end
        end
                   
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
        
        function clear(obj,src,~)
            if isfield(obj.handels,'traceEstimateTabPlot')
                    obj.deletetry(obj.handels.traceEstimateTabPlot)
             end
             if isfield(obj.handels,'analyzedTargets')
                for i=1:length(obj.handels.analyzedTargets)
                    obj.deletetry(obj.handels.analyzedTargets(i))
                end
             end
            if isfield(obj.handels,'spotSelectionRectangle')
                for i=1:length(obj.handels.spotSelectionRectangle)
                    obj.deletetry(obj.handels.spotSelectionRectangle{i})
                end
            end
            if isfield(obj.handels,'positionMarker') 
                obj.deletetry(obj.handels.positionMarker)
            end
            if isfield(obj.handels,'rastergramPlot') 
                obj.deletetry(obj.handels.rastergramPlot)
            end
            if isfield(obj.handels,'analyzedTargets') 
                for i=1:length(obj.handels.analyzedTargets)
                    obj.deletetry(obj.handels.analyzedTargets(i))
                end
            end
            if isfield(obj.handels,'imageStackSpotCam1') 
                obj.deletetry(obj.handels.imageStackSpotCam1)
            end
            if isfield(obj.handels,'imageStackSpotCam2') 
                obj.deletetry(obj.handels.imageStackSpotCam2)
            end

            obj.handels.selectedSpotIndex=[];
            if isfield(obj.handels,'traceEstimatePlot') 
                obj.deletetry(obj.handels.traceEstimatePlot)
            end
            if isfield(obj.handels,'imageStackCam1') 
                obj.deletetry(obj.handels.imageStackCam1)
            end
            if isfield(obj.handels,'imageStackCam2') 
                obj.deletetry(obj.handels.imageStackCam2)
            end
            if isfield(obj.handels,'timeLine') 
                obj.deletetry(obj.handels.timeLine)
            end
            if isfield(obj.handels,'initialClick') 
                obj.deletetry(obj.handels.initialClick)
            end
            if isfield(obj.handels,'spotLine') 
                obj.deletetry(obj.handels.spotLine)
            end
            
            if isfield(obj.handels,'componentPlot') 
                obj.deletetry(obj.handels.componentPlot)
            end
            if isfield(obj.handels,'componentAndRefPlot') 
                obj.deletetry(obj.handels.componentAndRefPlot)
            end
            if isfield(obj.handels,'overVieuwStack1') 
                obj.deletetry(obj.handels.overVieuwStack1)
            end
            if isfield(obj.handels,'overVieuwStack2') 
                obj.deletetry(obj.handels.overVieuwStack2)
            end
    
%             names = fieldnames(obj.handels)
%             for i=1:length(names)
%                 obj.deletetry(names{i});
%             end
            
            %results
            obj.handels.selectedSpotIndex=[];
            obj.addedvars = false;
            obj.AddedTrackingInfo = [];
            obj.sortMap=[];
            obj.spotNrs=[];
            obj.meanTime2FirstEvent=[];
            obj.rastergramStartframe =[];
            obj.meanDwellTime=[];
            obj.rawFitResultsCam1=[];
            obj.rawFitResultsCam2=[];
            obj.coordsCam1=[];
            obj.coordsCam2=[];
            obj.coord1=[];
            obj.coord2=[];
            obj.spotsIncluded=[];
            obj.maskFilt2=[];
            obj.maskFilt1=[];
            obj.photons1I=[];
            obj.photons2I=[];
            obj.bg2I=[];
            obj.bg1I=[];
            obj.deltax=[];
            obj.deltay=[];
            obj.seq=[];
            obj.Sigma1I=[];
            obj.intergrationTime=[];
            obj.status=0;
        end
        
        function changed = thresholdChange(obj,src,~)
            set(obj.handels.intensityTHSlider,'Value',obj.handels.Figure.UserData.sliderEventClass.TH(1));
            set(obj.handels.positionTHSlider,'Value',obj.handels.Figure.UserData.sliderEventClass.TH(2));
            set(obj.handels.backgroundRatioTHSlider,'Value',obj.handels.Figure.UserData.sliderEventClass.TH(3));
            set(obj.handels.sigmaTHSlider,'Value',obj.handels.Figure.UserData.sliderEventClass.TH(4));

            THPtemp = (obj.THPmax-1)/get(obj.handels.intensityTHSlider,'Max')*get(obj.handels.intensityTHSlider,'Value');
            THXtemp = (obj.THXmax-1)/get(obj.handels.positionTHSlider,'Max')*get(obj.handels.positionTHSlider,'Value');
            THBGtemp = (obj.THBGmax-1)/get(obj.handels.backgroundRatioTHSlider,'Max')*get(obj.handels.backgroundRatioTHSlider,'Value');
            THStemp = (obj.THSmax-1)/get(obj.handels.sigmaTHSlider,'Max')*get(obj.handels.sigmaTHSlider,'Value');
            if isempty(obj.THP) || isempty(obj.THX) || isempty(obj.THBG) || isempty(obj.THS) || (obj.THP ~= THPtemp || obj.THX ~= THXtemp || obj.THBG ~= THBGtemp || obj.THS ~= THStemp)              
                obj.THS = THStemp;
                obj.THP = THPtemp;
                obj.THX = THXtemp;
                obj.THBG = THBGtemp;
                if ~isempty(obj.NFrames)
                    obj.plotImages;
                end
                
                changed = true;
            else
                changed = false;
            end
        end
        
        function plotImages(obj,src,~)
            reDraw = ~isfield(obj.handels,'selectedSpotIndex') || isempty(obj.handels.selectedSpotIndex);

            if reDraw
                obj.handels.selectedSpotIndex =  get(obj.handels.spotSelectionSlider,'Value');
                maxSpots =  get(obj.handels.spotSelectionSlider,'Max');
                spotIdx = obj.handels.selectedSpotIndex:maxSpots:obj.NFrames*maxSpots;
                obj.handels.minimumSpotIntensity1 = min(dip_image(obj.rawFitResultsCam1.ROIStack(:,:,spotIdx)));
                obj.handels.maximumSpotIntensity1 = max(dip_image(obj.rawFitResultsCam1.ROIStack(:,:,spotIdx)));

                obj.handels.minimumSpotIntensity2 = min(dip_image(obj.rawFitResultsCam2.ROIStack(:,:,spotIdx)));
                obj.handels.maximumSpotIntensity2 = max(dip_image(obj.rawFitResultsCam2.ROIStack(:,:,spotIdx)));
            end
            
            set(obj.handels.spotIncludedCheckBox,'Value',obj.spotsIncluded(obj.handels.selectedSpotIndex));
            set(obj.handels.settingsTextTH,'String',sprintf('I=%0.2g, X=%0.2g, Ph x BG =%0.2g,S = %0.2g',obj.THP,obj.THX,obj.THBG, obj.THS));

            traces2_allfr = obj.getBinary(obj.handels.selectedSpotIndex);

            ymax = max(max(obj.bg1I(obj.handels.selectedSpotIndex,:,2-0)*(obj.PSFSigma*2)^2),max(max(obj.photons1I(obj.handels.selectedSpotIndex,:,2-0)),max(obj.photons1I(obj.handels.selectedSpotIndex,:,1+0))));
            if isfield(obj.handels,'traceEstimateTabPlot') && ~isempty(obj.handels.traceEstimateTabPlot)
                delete(obj.handels.traceEstimateTabPlot);
            end
            obj.handels.traceEstimateTabPlot = subplot(2,2,1:2,'Parent', obj.handels.traceEstimateTab,'buttondownfcn',@obj.clickline,'xlimmode','manual');
            
            axis([0 size(obj.photons1I,2) 0 ymax])

            line(1:size(obj.photons1I,2),obj.photons1I(obj.handels.selectedSpotIndex,:,2-0),'Color',obj.colorComplex,'Marker','x','parent',obj.handels.traceEstimateTabPlot)
            hold on
            if obj.tagetTrace
                line(1:size(obj.photons1I,2),obj.photons1I(obj.handels.selectedSpotIndex,:,1+0),'Color',obj.colorTarget,'parent',obj.handels.traceEstimateTabPlot,'tag','background','hittest','off')
            end
            line(1:size(obj.photons1I,2),obj.bg1I(obj.handels.selectedSpotIndex,:,2-0)*(obj.PSFSigma*2)^2,'Color',[.7 .5 0],'parent',obj.handels.traceEstimateTabPlot,'tag','background','hittest','off')
            line(1:size(obj.photons1I,2),traces2_allfr*ymax,'Color','b','parent',obj.handels.traceEstimateTabPlot,'tag','background','hittest','off')
            
            obj.parentFigure = obj.getParentFigure(obj.parent); 
            c = uicontextmenu('Parent', obj.parentFigure );
            obj.handels.traceEstimateTabPlot.UIContextMenu=c;
             
            uimenu(c,'Label','Export Figure','Callback',@obj.exportTraceFigure);
            uimenu(c,'Label','Export Data','Callback',@obj.exportTraceData);

            
            if ~reDraw & isfield(obj.handels,'timeLine') && ishandle(obj.handels.timeLine)
               xcoords = get(obj.handels.timeLine,'xdata');
               ii= xcoords(1);
            else            
                ii=obj.startFrame+10;
            end
            
            obj.handels.timeLine = line([ii ii],[0 ymax],'LineStyle','--','Color','k','Linewidth',1,'parent',obj.handels.traceEstimateTabPlot,'tag','timeLine','hittest','off');
            ylabel('I [# Photons]')
            xlabel(['Time x ' num2str(1./obj.intergrationTime) ' [s]'])
            drawnow

            if obj.tagetTrace
                legend('Component','Reference','Background','Binary','Orientation','horizontal') %,'Location','southoutside','Orientation','horizontal'
            else
                legend('Component','Background','Binary','Orientation','horizontal') %,'Location','southoutside','Orientation','horizontal'
            end

            hold off
         
            if (~isfield(obj.handels,'imageStackCam1') || isempty(obj.handels.imageStackCam1) || reDraw) && ( isempty(obj.draw3Dmovie) || obj.draw3Dmovie )
                x1=0.05;
                y1=0.1;

                obj.handels.imageStackCam1=subplot(1,2,2,'Parent', obj.handels.selectedSpotTab,'buttondownfcn',@(src,eventdata)obj.selectMarker(src,eventdata,1),'xlimmode','manual')

                params=[];

                obj.handels.overVieuwStack1 = imtool3D(obj.handels.imageStackCam1,@obj.showSliceCam1,params);
                set(obj.handels.overVieuwStack1.himg,'buttondownfcn',@(src,eventdata)obj.selectMarker(src,eventdata,1)) ;

                obj.handels.imageStackCam2=subplot(1,2,1,'Parent', obj.handels.selectedSpotTab,'buttondownfcn',@(src,eventdata)obj.selectMarker(src,eventdata,2),'xlimmode','manual');

                obj.handels.overVieuwStack2 = imtool3D(obj.handels.imageStackCam2,@obj.showSliceCam2,params);
                set(obj.handels.overVieuwStack2.himg,'buttondownfcn',@(src,eventdata)obj.selectMarker(src,eventdata,2));

%             else
%                 obj.setSelRecPositions;
            end
            obj.plotTimePoint;  
        end    
        
        function showTargetTrace(obj,bol)
           obj.tagetTrace=bol;
        end       
        
        function setSelRecPositions(obj,index)
            % parent object should be a im3d object!
            
            if nargin < 2
              index = [1,2];
            end
          
            boxsize=10;
            x={{},{}}; y={{},{}};
            if any(index ==1) && ~isempty(obj.coord1)
                x{1}=obj.coord1(obj.handels.selectedSpotIndex,1);
                y{1}=obj.coord1(obj.handels.selectedSpotIndex,2);
            end
            if any(index ==2) && ~isempty(obj.coord2)
                x{2}=obj.coord2(obj.handels.selectedSpotIndex,1);
                y{2}=obj.coord2(obj.handels.selectedSpotIndex,2);
            end
             
            for i=index
                if ~isempty(x{i})
                    if isfield(obj.handels,'spotSelectionRectangle') && ~isempty(obj.handels.spotSelectionRectangle) && i <= size(obj.handels.spotSelectionRectangle,2) && ~isempty(obj.handels.spotSelectionRectangle{i}) && ishandle(obj.handels.spotSelectionRectangle{i})
                        set(obj.handels.spotSelectionRectangle{i},'Position',[x{i}-boxsize/2+1 y{i}-boxsize/2+1 boxsize boxsize]);
                     else 
                         obj.handels.spotSelectionRectangle{i} = rectangle('Position',[x{i}-boxsize/2+1 y{i}-boxsize/2+1 boxsize boxsize],'lineWidth',2,'EdgeColor','y','tag','selectionBox','hittest','off','buttondownfcn',@(src,eventdata)obj.selectMarker(src,eventdata,i))
                    end
                end
            end
        end
        
        function [img,  stackLength] = showSliceCam1(obj,im3d,handle,newSlice,update)
            if nargin <= 4
                update = true;
            end
            if nargin < 3
                newSlice=obj.startFrame;
                handle=[];
            end
            title(sprintf('t=%d',round(newSlice-obj.startFrame)),'parent',obj.handels.imageStackCam1)

            
            stackLength = max(1,size(obj.dataObject.getData(1),3));
            obj.colorTarget;
            img = repmat(obj.dataObject.getData(1,[],[],newSlice),[1 1 3]);
          
            if ~isempty(handle) && ishandle(handle) && update
                switch im3d.imgIntensityMapping
%                     case 'linear'
%                       handle.CData = double(stretch(()));    
%                     case 'log'
%                         handle.CData = double(stretch(real(log(1+obj.dataObject.dataCam{1}(:,:,newSlice)))));
                        
                        case 'linear'
                            a = sort(reshape(obj.dataObject.getData(1,[],[],newSlice),[1 size(obj.dataObject.getData(1),1)*size(obj.dataObject.getData(1),1)]));
                             imgdata = double(imadjust(mat2gray((obj.dataObject.getData(1,[],[],newSlice))),[  max(0,a(max(1,ceil(size(a,2)*im3d.pStretch.low)))./max(a)) min(1,a(min(end,round(size(a,2)*im3d.pStretch.high)))./max(a))],[0 1]));
                             handle.CData =  cat(3,obj.colorTarget(1).*imgdata,obj.colorTarget(2).*imgdata,obj.colorTarget(3).*imgdata);
                        case 'log'
                            imglog = (real(log(1+obj.dataObject.getData(1,[],[],newSlice))));
                            a = sort(reshape(imglog,[1 size(obj.dataObject.getData(1),1)*size(obj.dataObject.getData(1),1)]));
                            imgdata = double(stretch(imglog,max(0,[  im3d.pStretch.low]),[ min(1,im3d.pStretch.high)*100],[0],[1]));
                            handle.CData =  cat(3,obj.colorTarget(1).*imgdata,obj.colorTarget(2).*imgdata,obj.colorTarget(3).*imgdata);
                end
                colormap(handle.Parent,'gray')    
                idx = obj.coordsCam1(:,3)+1 == newSlice; %max(obj.startFrame,newSlice);
                if sum(idx) ~= 0
                    obj.coord1(:,1) = obj.coordsCam1(idx,1);
                    obj.coord1(:,2) = obj.coordsCam1(idx,2);
                    obj.updateBoxes;
                    setSelRecPositions(obj,1);
                 else
                    obj.coord1(:,1) = obj.segmentedSpots(:,1);
                    obj.coord1(:,2) = obj.segmentedSpots(:,2);
                    obj.updateBoxes;
                    setSelRecPositions(obj,1);                                
                end
            end
        end

        function [img,  stackLength] = showSliceCam2(obj,im3d,handle,newSlice,update)
             if nargin <= 4
                update = true;
            end
            if nargin < 3
                newSlice=obj.startFrame;
                handle=[];
            end
            title(sprintf('t=%d',round(newSlice-obj.startFrame)),'parent',obj.handels.imageStackCam2)


            stackLength = max(1,size(obj.dataObject.getData(2),3));
            img = repmat(obj.dataObject.getData(2,[],[],newSlice),[1 1 3]);
            

            if ~isempty(handle) && ishandle(handle) && update
                switch im3d.imgIntensityMapping
                         case 'linear'
                            a = sort(reshape(obj.dataObject.getData(2,[],[],newSlice),[1 size(obj.dataObject.getData(2),1)*size(obj.dataObject.getData(2),1)]));
                            imgdata= double(imadjust(mat2gray((obj.dataObject.getData(2,[],[],newSlice))),[  max(0,a(max(1,ceil(size(a,2)*im3d.pStretch.low)))./max(a)) min(1,a(min(end,round(size(a,2)*im3d.pStretch.high)))./max(a))],[0 1]));
                            
                            handle.CData =  cat(3,obj.colorComplex(1).*imgdata,obj.colorComplex(2).*imgdata,obj.colorComplex(3).*imgdata);
                        case 'log'
                            imglog = (real(log(1+obj.dataObject.getData(2,[],[],newSlice))));
                            a = sort(reshape(imglog,[1 size(obj.dataObject.getData(2),1)*size(obj.dataObject.getData(2),1)]));
                            imgdata = double(stretch(imglog,[  max(0,im3d.pStretch.low)],[ min(100,im3d.pStretch.high*100)],[0],[1]));
                            handle.CData =  cat(3,obj.colorComplex(1).*imgdata,obj.colorComplex(2).*imgdata,obj.colorComplex(3).*imgdata);
                end
            colormap(handle.Parent,'gray')
            idx = obj.coordsCam1(:,3)+1 == newSlice; %max(obj.startFrame,newSlice);
            if sum(idx) ~= 0
                obj.coord2(:,1) = obj.coordsCam2(idx,1);
                obj.coord2(:,2) = obj.coordsCam2(idx,2);
                obj.updateBoxes;
                setSelRecPositions(obj,2);
            else
                icoordsCam2=transformPointsInverse(obj.alignmentObject.tformTotal,obj.segmentedSpots); 
                obj.coord2(:,1) = icoordsCam2(:,1);
                obj.coord2(:,2) = icoordsCam2(:,2);
                obj.updateBoxes;
                setSelRecPositions(obj,2);                
            end
            end
        end

        function selectMarker(obj,src,ev,index)
            clicked=get(gca,'currentpoint');
            xcoord=clicked(1,1,1);
            ycoord=clicked(1,2,1);
            if nargin < 4 || isempty(index) || index ==1
                obj.handels.selectedSpotIndex = knnsearch(obj.coord1,[xcoord ycoord]);
            elseif index == 2
                obj.handels.selectedSpotIndex = knnsearch(obj.coord2,[xcoord ycoord]);
            end
            
            switch get(gcf,'selectiontype')
                case 'normal'
                    setSelRecPositions(obj);
                    obj.plotImages;
                case {'open','extend'}
                    if ~isempty(obj.spotsIncluded) && ~isempty(obj.handels.selectedSpotIndex)
                        obj.spotsIncluded(obj.handels.selectedSpotIndex) = ~obj.spotsIncluded(obj.handels.selectedSpotIndex);
                    end
                    obj.updateBoxes;
            end
          end
        function fig = getParentFigure(~,fig)
            while ~isempty(fig) && ~strcmp('figure', get(fig,'type'))
              fig = get(fig,'parent');
            end
        end
        function clickline(obj,src,ev)
           if strcmp(get(gcf,'Selectiontype'),'alt')
                a = get(gcf,'currentpoint');
                obj.handels.traceEstimateTabPlot.UIContextMenu.Position = a(1,1:2);
                obj.handels.traceEstimateTabPlot.UIContextMenu.Visible='on';
            else
                clicked=get(gca,'currentpoint');
                obj.handels.initialClick=clicked(1,1,1);

                set(gcf,'windowbuttonmotionfcn',@obj.dragline);
                set(gcf,'windowbuttonupfcn',@obj.dragdone);
           end
        end
        
        function dragline(obj,src,ev)
          
            clicked=get(gca,'currentpoint');
            xcoord=clicked(1,1,1);
            set(obj.handels.timeLine,'xdata',[xcoord xcoord]);
            ii = round(xcoord);
            if ~isempty(obj.seq) && ii <= size(obj.seq,2) && ii>=1
                title(sprintf('spot ID = %d\n Frame = %d, \\Delta X = %0.2g, I = %0.2g, bg=%0.2g, S = %0.2g',obj.handels.selectedSpotIndex,ii,obj.seq(obj.handels.selectedSpotIndex,ii),obj.photons1I(obj.handels.selectedSpotIndex,ii,2-0),obj.bg1I(obj.handels.selectedSpotIndex,ii,2-0)*(obj.PSFSigma*2)^2,obj.Sigma1I(obj.handels.selectedSpotIndex,ii,2-0)),'Parent',obj.handels.traceEstimateTabPlot);
            end            
        end

        function dragdone(obj,src,ev)
            clicked=get(gca,'currentpoint');
            finalClick=clicked(1,1,1);
            
            if finalClick == obj.handels.initialClick
                xcoord = get(obj.handels.timeLine,'xdata')
                if finalClick < xcoord(1)
                    set(obj.handels.timeLine,'xdata',xcoord-1);
                else
                    set(obj.handels.timeLine,'xdata',xcoord+1);
                end
            end

            set(gcf,'windowbuttonmotionfcn','');
            set(gcf,'windowbuttonupfcn','');
            notify(obj,'timePointChanged');
        end
              
        function [meanDwellTime, meanTime2FirstEvent,cia] = calculateRastergram(obj,traces2_allfr,ttb,allSpotsOverwrite,sortState)
            obj.classAdded
            if nargin < 5 || isempty(sortState)
                sortState = 1;
            end
            
            if ~obj.status
                error('Analysis has to be completed before rastergram can be calculated!')
            end
            
            if isempty(obj.numberOfSpotsForRastergram) || (nargin > 3 && ~isempty(allSpotsOverwrite) && allSpotsOverwrite)
                idx = true(size(traces2_allfr,1),1);
            else
                idx = false(size(traces2_allfr,1),1);
                idx(randperm(size(traces2_allfr,1),min(size(traces2_allfr,1),obj.numberOfSpotsForRastergram))) = true;
            end
            obj.spotNrs = find(obj.spotsIncluded);
            obj.spotNrs = obj.spotNrs(idx);
            traces2_allfr= traces2_allfr(idx,:);

            if isfield(obj.handels,'rastergramPlot')
                obj.deletetry(obj.handels.rastergramPlot);
            end
            obj.handels.rastergramPlot = subplot(1,1,1,'Parent', obj.handels.rastergramTap);
            obj.handels.rastergramPlot.ButtonDownFcn = @obj.clickhline;
            obj.handels.rastergramPlot.XLimMode = 'manual';
            
            [obj.rasterGramOuput.evi,obj.rasterGramOuput.cev,obj.rasterGramOuput.cia,obj.meanTime2FirstEvent,obj.meanDwellTime,obj.handels.rastergramImg] = newRastergram(obj.handels.rastergramPlot,traces2_allfr,ttb,obj.intergrationTime,obj.align,obj.delete,obj.sortType,sortState);
            obj.handels.rastergramPlot.UserData.plotData = obj.rasterGramOuput.cev;
            
            
            obj.parentFigure = obj.getParentFigure(obj.parent); 
            c = uicontextmenu('Parent', obj.parentFigure );
            obj.handels.rastergramPlot.UIContextMenu=c;
             
            uimenu(c,'Label','Export Figure','Callback',@obj.exportRastergramFig);
            uimenu(c,'Label','Export Data','Callback',@obj.exportRastergramData);
            
            obj.sortMap=obj.rasterGramOuput.evi;
            obj.spotNrMap=obj.spotNrs(obj.sortMap)
            obj.handels.rastergramImg.ButtonDownFcn = @obj.clickhline;
            if sum(obj.handels.selectedSpotIndex == obj.spotNrs) == 0
                 obj.handels.selectedSpotIndex = obj.spotNrs(round(size(obj.spotNrs,2)/2));
                 obj.plotImages;
            end
            hold on

            idx = obj.handels.selectedSpotIndex == obj.spotNrMap;
            if ~isempty(idx)
                obj.handels.spotLine = line(obj.handels.rastergramImg.Parent.XLim,[1 1].*find(idx),'LineStyle','--','Color','k','Linewidth',2,'parent',obj.handels.rastergramImg.Parent,'tag','timeLine');
                obj.handels.spotLine.ButtonDownFcn=@obj.clickhline;
            end
            title(sprintf('spot ID = %d, meanTime2FirstEvent = %0.2g, meanDwellTime = %0.2g',obj.handels.selectedSpotIndex,obj.meanTime2FirstEvent,obj.meanDwellTime),'Parent', obj.handels.rastergramPlot)
            meanTime2FirstEvent = obj.meanTime2FirstEvent;
            meanDwellTime = obj.meanDwellTime;
        end
        
        function [traces2_allfr,data] = getBinary(obj,spotsIncluded)
            if ~isempty(spotsIncluded)
               traces2_allfr = obj.photons1I(spotsIncluded,:,2-0) > obj.THP...
                & sqrt(obj.deltay(spotsIncluded,:).^2+obj.deltax(spotsIncluded,:).^2) < obj.THX...
                & obj.photons1I(spotsIncluded,:,2-0) > obj.THBG*obj.bg1I(spotsIncluded,:,2-0)*(obj.PSFSigma*2)^2 ...
                & obj.Sigma1I(spotsIncluded,:,2-0) < obj.THS;  
                if nargout > 1
                    data=cat(3,obj.deltay(spotsIncluded,:).^2+obj.deltax(spotsIncluded,:).^2,obj.Sigma1I(spotsIncluded,:,2-0),obj.bg1I(spotsIncluded,:,2-0)*(obj.PSFSigma*2)^2,obj.photons1I(spotsIncluded,:,2-0));
                end
            else
                traces2_allfr=[];
                data=[];
            end
        end
        
        function [traces2_allfr, ttb] = getTraceMatrix(obj)
            traces2_in = obj.getBinary(obj.spotsIncluded);

            traces2_allfr=[];
            if isfield(obj.AddedTrackingInfo,'spots')
                idx = [0 cumsum(obj.AddedTrackingInfo.spots)];
            else
                idx=[0 size(traces2_in,1)];
            end
            
            for i=1:length(obj.rastergramStartframe)
                dT = max(obj.rastergramStartframe) - obj.rastergramStartframe(i);
                traces2_allfr = cat(1,traces2_allfr,traces2_in(1+idx(i):idx(i+1),obj.rastergramStartframe(i):size(traces2_in,2)-dT));
            end
            
            obj.spotNrs = find(obj.spotsIncluded);
           
            if isempty(obj.intergrationTime)
                obj.intergrationTime=1;
            end
            if ischar(obj.intergrationTime)
                ws = load(obj.intergrationTime,'vid')
                ttb = ws.vid.ttb;
            elseif isnumeric(obj.intergrationTime)
                ttb=1000*[0:obj.intergrationTime:obj.intergrationTime*size(traces2_allfr,2)]; % This is an array containing relative time stamps for all frames
            end
        end
        
        function time21st = getTimes2FirstEvent(obj,calculateRastergram,dataSetIndex)
            if nargin < 2 || isempty(calculateRastergram)
                calculateRastergram = true;
            end
            
            if calculateRastergram
                [traces2_allfr, ttb] = obj.getTraceMatrix;
                if nargin > 2 && ~isempty(dataSetIndex)
                    idx = [0 cumsum(obj.AddedTrackingInfo.spots)];
                    obj.calculateRastergram(traces2_allfr(1+idx(dataSetIndex):idx(dataSetIndex+1),:),ttb);
                else
                    obj.calculateRastergram(traces2_allfr,ttb);
                end
            end
            
            rastStruct = obj.rasterGramOuput;
            logik=[];
            logik=(rastStruct.cia(:,1)==-2);
            Time2FirstEvent=rastStruct.cia(logik,5); 
            time_empty=max(Time2FirstEvent)*ones(length(Time2FirstEvent),1);
            t=(Time2FirstEvent~=time_empty);

            time21st=Time2FirstEvent(t);
        end
        
        function dwellTimes = getDwellTimes(obj,calculateRastergram,dataSetIndex)
            if nargin < 2 || isempty(calculateRastergram)
                calculateRastergram = true;
            end
            
            if calculateRastergram
                [traces2_allfr, ttb] = obj.getTraceMatrix;
                if nargin > 2 && ~isempty(dataSetIndex)
                    idx = [0 cumsum(obj.AddedTrackingInfo.spots)];
                    obj.calculateRastergram(traces2_allfr(1+idx(dataSetIndex):idx(dataSetIndex+1),:),ttb);
                else
                    obj.calculateRastergram(traces2_allfr(:,:),ttb);
                end
            end

            rastStruct = obj.rasterGramOuput;
            logik=[];
            logik=rastStruct.cia(:,1)==1;
            dwellTimes=rastStruct.cia(logik,5);
        end
     
        
        function setRastergramSettings(obj,src,~)
             Title = 'Rastergram Settings';

            %%%% SETTING DIALOG OPTIONS
            % Options.WindowStyle = 'modal';
            Options.Resize = 'on';
            Options.Interpreter = 'tex';
            Options.CancelButton = 'on';
            Options.ApplyButton = 'off';
            % Options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration
            Option.Dim = 1; % Horizontal dimension in fields

            Prompt = {};
            Formats = {};
            % DefAns = struct([]);


            Prompt(1,:) = {'Startframe Rastergrams','rastergramStartframe',[]};
            Formats(1,1).type = 'edit';
            Formats(1,1).format = 'text';
            Formats(1,1).size = [-1 0];
            Formats(1,1).span = [1 1];  % item is 1 field x 3 fields
            if ~isscalar(obj.rastergramStartframe)
                DefAns.rastergramStartframe = [ '[' num2str(obj.rastergramStartframe) ']' ]; %[pwd '/Dark.tif'];
            else
                 DefAns.rastergramStartframe = num2str(obj.rastergramStartframe); %[pwd '/Dark.tif'];
            end

            Prompt(2,:) = {'Sorting','List',[]};
            Formats(2,1).type = 'list';
            Formats(2,1).style = 'popupmenu';
            Formats(2,1).items = {'none','duration of first','duration of last','binding of first'};
            
            Prompt(3,:) = {'Align' 'align',[]};
            Formats(3,1).type = 'list';
            Formats(3,1).style = 'popupmenu';
            Formats(3,1).items = {'none','first binding event','last departure event'};
                        
            Prompt(4,:) = {'Delete eventless' 'delete',[]};
            Formats(4,1).type = 'check';
            DefAns.delete = obj.delete;

            Prompt(5,:) = {'Number of Spots','numberOfSpotsForRastergram',[]};
            Formats(5,1).type = 'edit';
            Formats(5,1).format = 'text';
            Formats(5,1).size = [-1 0];
            Formats(5,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.numberOfSpotsForRastergram = num2str(obj.numberOfSpotsForRastergram); 

            [result,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);

            if ~Cancelled
                obj.sortType = result.List-1;
                obj.align = result.align-1;
                obj.delete = result.delete;
                obj.rastergramStartframe = eval(result.rastergramStartframe);
                if isempty(result.numberOfSpotsForRastergram)
                    obj.numberOfSpotsForRastergram =[];
                else
                    obj.numberOfSpotsForRastergram = str2num(result.numberOfSpotsForRastergram);
                end
            end
        end
        
        function bootStrapRastergramMenu(obj,src,~)
              if ~isempty(obj.coordsCam1)
                prompt={'Enter number of bootstraps'};
                name = 'Enter number of bootstraps';

                defaultans = {num2str(10)};
                options.Interpreter = 'tex';

                answer = inputdlg(prompt,name,[1 40],defaultans,options);
                if ~isempty(answer)
                    obj.Nbootstraps = str2num(answer{1});  

                    [meanDwellTimeBS,meanTime2FirstEventBS] = obj.bootStrapRastergram(obj.Nbootstraps);
                end
              else
                  h = msgbox('Perform analysis first')
              end
        end
        
        function [meanDwellTimeBS,meanTime2FirstEventBS] = bootStrapRastergram(obj,N,sortState)
            
            if nargin < 3 || isempty(sortState)
                sortState = 1;
            end
            
            if ~obj.status
                error('Analysis has to be completed before rastergram can be calculated!')
            end
            
            if nargin < 2 || isempty(N)
                N=10;
            end
            obj.Nbootstraps = N;
            
            [traces2_allfr, ttb] = getTraceMatrix(obj);

            meanDwellTimeBS=zeros(N,1);
            meanTime2FirstEventBS=zeros(N,1);
            [~, ~,obj.rasterGramOuput.cia] = calculateRastergram(obj,traces2_allfr,ttb,true,sortState);
            
            for i=1:N
                idx = randperm(size(traces2_allfr,1),round(size(traces2_allfr,1)*0.9));
                
                sel = zeros(size(obj.rasterGramOuput.cia,1),1);
                sel(idx)=1;
            
                logikTime21st=(obj.rasterGramOuput.cia(:,1)==-2) & sel;
                if any(logikTime21st)
                    Time2FirstEvent=obj.rasterGramOuput.cia(logikTime21st,5); 
                    time_empty=max(Time2FirstEvent)*ones(length(Time2FirstEvent),1);
                    t=(Time2FirstEvent~=time_empty);
                    time21st=Time2FirstEvent(t);
                    meanTime2FirstEventBS(i) = mean(time21st);
                end
                
                logikDwellTime=(obj.rasterGramOuput.cia(:,1)==1) & sel;
                if any(logikDwellTime)
                    DwellTime=obj.rasterGramOuput.cia(logikDwellTime,5);    
                    meanDwellTimeBS(i) = mean(DwellTime);
                end
                
                title(sprintf('spot ID = %d, meanTime2FirstEvent = %0.2g \\pm %0.2g, meanDwellTime = %0.2g \\pm%0.2g, At %d of %d',obj.handels.selectedSpotIndex,mean(meanTime2FirstEventBS),...
                std(meanTime2FirstEventBS),mean(meanDwellTimeBS),std(meanDwellTimeBS),i,N),'Parent', obj.handels.rastergramPlot)
                
            drawnow;
            end
        end
        
        function runRastergram(obj,src,~)
             if ~obj.status
                error('Analysis has to be completed before rastergram can be calculated!')
            end
            [traces2_allfr, ttb] = getTraceMatrix(obj);
            if ~isempty(traces2_allfr)
                obj.calculateRastergram(traces2_allfr,ttb);
            else
                h = msgbox('No traces', 'Error','error');
            end
        end
        
        
        
          function plotTimePoint(obj,src,~)

            xcoords = get(obj.handels.timeLine,'xdata');
            ii = round(xcoords(1));

            frames2 = obj.rawFitResultsCam2.Frame(obj.maskFilt2,:);
            frames1 = obj.rawFitResultsCam1.Frame(obj.maskFilt1,:);
    
            idx = find(frames2==ii+min(frames2)-1);
            A(:,:,1+obj.flipCameras) = squeeze(obj.rawFitResultsCam2.ROIStack(:,:,idx(obj.handels.selectedSpotIndex)));
      
            idx1 = find(frames1==ii+min(frames2)-1);
            A(:,:,2-obj.flipCameras) = squeeze(obj.rawFitResultsCam1.ROIStack(:,:,idx1(obj.handels.selectedSpotIndex)));
            
            xxyy =  squeeze(obj.coordsComplex(obj.handels.selectedSpotIndex,ii,:));
            xxyytarget =  squeeze(obj.coordsTarget(obj.handels.selectedSpotIndex,ii,:));
            
            if ~isempty(obj.seq)
                b = obj.seq(obj.handels.selectedSpotIndex,:);
                title(sprintf('spot ID = %d\n Frame = %d, \\Delta X = %0.2g, I = %0.2g, bg=%0.2g, S=%0.2g',obj.handels.selectedSpotIndex,ii,sqrt(sum((xxyy-xxyytarget).^2)),obj.photons1I(obj.handels.selectedSpotIndex,ii,2-0),obj.bg1I(obj.handels.selectedSpotIndex,ii,2-0)*(obj.PSFSigma*2)^2,obj.Sigma1I(obj.handels.selectedSpotIndex,ii,2-0)),'Parent',obj.handels.traceEstimateTabPlot);
            end              


            gimr = mat2gray(A(:,:,1)',[min(obj.handels.minimumSpotIntensity1,obj.handels.minimumSpotIntensity2) max(obj.handels.maximumSpotIntensity1,obj.handels.maximumSpotIntensity2)]); %flip(flip(A(:,:,ii+1),1),2));
            gimr = uint8(256.*cat(3,obj.colorComplex(1).*gimr, obj.colorComplex(2).*gimr, obj.colorComplex(3).*gimr)); %flip(flip(A(:,:,ii+1),1),2));      
            
            if isfield(obj.handels,'imageStackSpotCam1') && ~isempty(obj.handels.imageStackSpotCam1) && isvalid(obj.handels.imageStackSpotCam1)
                obj.handels.imageStackSpotCam1.CData=gimr;
            else
                obj.handels.componentPlot = subplot(2,2,3,'Parent', obj.handels.traceEstimateTab);
                obj.handels.imageStackSpotCam1 = subimage(gimr);
                axis off;
                hold on
                title('Component','parent',obj.handels.componentPlot)
            end

            gimg = mat2gray(flip(flip(A(:,:,2),1),2),[min(obj.handels.minimumSpotIntensity1,obj.handels.minimumSpotIntensity2) max(obj.handels.maximumSpotIntensity1,obj.handels.maximumSpotIntensity2)]); 
            gimg = uint8(gimg * 256);
            gimg = cat(3,obj.colorTarget(1).*gimg,obj.colorTarget(2).*gimg,obj.colorTarget(3).*gimg);

            rgbImage = gimr+gimg;
            
            if isfield(obj.handels,'componentAndRefPlot') && ~isempty(obj.handels.componentAndRefPlot) && isvalid(obj.handels.componentAndRefPlot)
                obj.handels.imageStackSpotCam2.CData=rgbImage;
                 obj.handels.imageStackSpotCam2.CDataMapping='scaled';
            else
                obj.handels.componentAndRefPlot = subplot(2,2,4,'Parent', obj.handels.traceEstimateTab);
                obj.handels.imageStackSpotCam2 = subimage(rgbImage);
                axis off
                xlim([1 size(rgbImage,1)])
                ylim([1 size(rgbImage,2)])
                title('Reference & Component','parent',obj.handels.componentAndRefPlot)
            end
            drawnow
            hold on
            
            if isfield(obj.handels, 'positionMarker') && length(obj.handels.positionMarker) >= 1 && ishandle(obj.handels.positionMarker(1))
                delete(obj.handels.positionMarker(1))
            end
            obj.handels.positionMarker(1) = plot(-0.1,-0.1,'h','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,'parent',obj.handels.componentAndRefPlot)
           
            
            if isfield(obj.handels,'positionMarker') && length(obj.handels.positionMarker) >= 2 && ishandle(obj.handels.positionMarker(2))
                delete(obj.handels.positionMarker(2))
            end
            obj.handels.positionMarker(2) = plot(xxyy(1)+1,xxyy(2)+1,'h','MarkerFaceColor',obj.colorComplex,'MarkerEdgeColor','k','MarkerSize',8,'parent',obj.handels.componentAndRefPlot)
          
              if obj.colocalized
                  if length(obj.handels.positionMarker) >= 3 && ishandle(obj.handels.positionMarker(3))
                        delete(obj.handels.positionMarker(3))
                  end
                  obj.handels.positionMarker(3) = plot(xxyytarget(1)+1,xxyytarget(2)+1,'h','MarkerFaceColor',obj.colorTarget,'MarkerEdgeColor','k','MarkerSize',8,'parent',obj.handels.componentAndRefPlot)
              else
                  if length(obj.handels.positionMarker) >= 3 && ishandle(obj.handels.positionMarker(3))
                        delete(obj.handels.positionMarker(3))
                  end
                   nocoloxxyy = -obj.rawFitResultsCam2.RoiStart(idx(obj.handels.selectedSpotIndex),:)+obj.coord1((obj.handels.selectedSpotIndex),:);
                   obj.handels.positionMarker(3) = plot(xxyytarget(1)+1,xxyytarget(2)+1,'h','MarkerFaceColor',obj.colorTarget,'MarkerEdgeColor','k','MarkerSize',8,'parent',obj.handels.componentAndRefPlot)
              end

                if length(obj.handels.positionMarker) >= 4 &&ishandle(obj.handels.positionMarker(4))
                    delete(obj.handels.positionMarker(4))
                end
            obj.handels.positionMarker(4) = plot(xxyy(1)+1,xxyy(2)+1,'h','MarkerFaceColor',obj.colorComplex,'MarkerEdgeColor','k','MarkerSize',8,'parent',obj.handels.componentPlot)
            if obj.colocalized
                if length(obj.handels.positionMarker) >= 5 &&ishandle(obj.handels.positionMarker(5))
                    delete(obj.handels.positionMarker(5))
                end
                obj.handels.positionMarker(5) = plot(xxyytarget(1)+1,xxyytarget(2)+1,'h','MarkerFaceColor',obj.colorTarget,'MarkerEdgeColor','k','MarkerSize',8,'parent',obj.handels.componentPlot)
            else
                 if length(obj.handels.positionMarker) >= 5 &&ishandle(obj.handels.positionMarker(5))
                    delete(obj.handels.positionMarker(5))
                end
                
                nocoloxxyy = -obj.rawFitResultsCam2.RoiStart(idx(obj.handels.selectedSpotIndex),:)+obj.coord1((obj.handels.selectedSpotIndex),:);
                obj.handels.positionMarker(5) = plot(xxyytarget(1)+1,xxyytarget(2)+1,'h','MarkerFaceColor',obj.colorTarget,'MarkerEdgeColor','k','MarkerSize',8,'parent',obj.handels.componentPlot)

            end
          end
        
        function clickhline(obj,src,ev)               
            if strcmp(get(gcf,'Selectiontype'),'alt')
                a = get(gcf,'currentpoint');
                obj.handels.rastergramPlot.UIContextMenu.Position = a(1,1:2);
                obj.handels.rastergramPlot.UIContextMenu.Visible='on';
            else
                clicked=get(gca,'currentpoint');
                obj.handels.initialClick=clicked(1,2,1);
                set(gcf,'windowbuttonmotionfcn',@obj.draghline);
                set(gcf,'windowbuttonupfcn',@obj.dragdonehline);
            end
        end
        
        function draghline(obj,src,ev)
          
            clicked=get(gca,'currentpoint');
            ycoord=clicked(1,2,1);
            obj.handels.spotLine.YData = [ycoord ycoord];
            drawnow;
            if ~isempty(obj.seq) && any(obj.sortMap>0)
                 obj.handels.selectedSpotIndex = obj.spotNrMap(min(size(obj.spotNrMap,1),max(1,round(ycoord+0.5))));
                 title(sprintf('spot ID = %d, meanTime2FirstEvent = %0.2g, meanDwellTime = %0.2g',obj.handels.selectedSpotIndex,obj.meanTime2FirstEvent,obj.meanDwellTime),'Parent', obj.handels.rastergramPlot)
            end       
        end

        function dragdonehline(obj,src,ev)
            clicked=get(gca,'currentpoint');
            finalClick=clicked(1,2,1);
            
            if finalClick == obj.handels.initialClick
                ycoord = get(obj.handels.spotLine,'ydata');

                if finalClick < ycoord(1)
                    obj.handels.selectedSpotIndex = obj.spotNrMap(round(round(ycoord(1))-1));
                    set(obj.handels.spotLine,'ydata',ycoord-1);
                else
                    obj.handels.selectedSpotIndex = obj.spotNrMap(round(round(ycoord(1))+1));
                    set(obj.handels.spotLine,'ydata',ycoord+1);
                end
                title(sprintf('spot ID = %d, meanTime2FirstEvent = %0.2g, meanDwellTime = %0.2g',obj.handels.selectedSpotIndex,obj.meanTime2FirstEvent,obj.meanDwellTime),'Parent', obj.handels.rastergramPlot)
            end

            set(gcf,'windowbuttonmotionfcn','');
            set(gcf,'windowbuttonupfcn','');
            obj.plotImages;
        end
        
        function vs = savevars(obj,src,~)

            if ~obj.addedvars
                class = obj.classSingle;
            else
                class = obj.classAdded;
            end
            vs.class=class;
            vs.PSFSigma=obj.PSFSigma;
            vs.startFrame = obj.startFrame;
            vs.detectionFrame = obj.detectionFrame;
            vs.NFrames = obj.NFrames;
            vs.flipCameras = obj.flipCameras;
            vs.subRegionSelection = obj.subRegionSelection;
            vs.colocalized = obj.colocalized;
            vs.rastergramStartframe=obj.rastergramStartframe;
            vs.numberOfSpotsForRastergram=obj.numberOfSpotsForRastergram;
            vs.THP=obj.THP;
            vs.THX=obj.THX;
            vs.THBG=obj.THBG;
            vs.THS=obj.THS;
            vs.Nbootstraps=obj.Nbootstraps;
            vs.tagetTrace =obj.tagetTrace;
            vs.colorTarget = obj.colorTarget;
            vs.colorComplex= obj.colorComplex;

        
            %results
            vs.sortMap=obj.sortMap;
            vs.spotNrs=obj.spotNrs;
            vs.meanTime2FirstEvent=obj.meanTime2FirstEvent;
            vs.meanDwellTime=obj.meanDwellTime;
            vs.rawFitResultsCam1=obj.rawFitResultsCam1;
            vs.rawFitResultsCam2=obj.rawFitResultsCam2;
            vs.coordsCam1=obj.coordsCam1;
            vs.coordsCam2=obj.coordsCam2;
            vs.THPmax=obj.THPmax;
            vs.THXmax=obj.THXmax;
            vs.THSmax=obj.THSmax;
            vs.THBGmax=obj.THBGmax;
            
    
            vs.coord1=obj.coord1;
            vs.coord2=obj.coord2;
            vs.spotsIncluded=obj.spotsIncluded;
            vs.maskFilt2=obj.maskFilt2;
            vs.maskFilt1=obj.maskFilt1;
            vs.photons1I=obj.photons1I;
            vs.photons2I=obj.photons2I;
            vs.bg2I=obj.bg2I;
            vs.bg1I=obj.bg1I;
            vs.Sigma1I=obj.Sigma1I;
            vs.deltax=obj.deltax;
            vs.deltay=obj.deltay;
            vs.seq=obj.seq;
            vs.intergrationTime=obj.intergrationTime;
            vs.segmentedSpots=obj.segmentedSpots;
            vs.coordsComplex=obj.coordsComplex;
            vs.coordsTarget=obj.coordsTarget;

           if nargin > 1 && isobject(src)
                [file, pathname] = uiputfile([class '.mat'] );
                if file ~= 0
                    answer{1} = fullfile(pathname, file);
                else
                    answer=[];
                end
                if ~isempty(answer)
                    save(answer{1},'vs')
                end
            end
        end

        function ws = loadvars(obj,src,updateplot)
            
            obj.clear;
            
            if ~obj.addedvars
                class = obj.classSingle;
            else
                class = obj.classAdded;
            end
            if ishandle(src)
%                 prompt={'File name'};
%                 name = 'Load';
%                 defaultans = {['guiAnalyzeDynamics' date '.mat']};
%                 options.Interpreter = 'tex';
%                 answer = inputdlg(prompt,name,[1 40],defaultans,options);

                [file, pathname] = uigetfile([class date '.mat'] );
                if file ~= 0
                    answer{1} = fullfile(pathname, file);
                else
                    answer=[];
                end
            elseif ischar(src)
                answer{1} = src;
            else
                answer{1} = [];  
            end
            if ~isempty(answer)
                
                obj.AddedTrackingInfo.files=[];
                obj.AddedTrackingInfo.files{1} = answer;

                if nargin < 3 || isempty(updateplot) || isobject(updateplot)
                     updateplot=true;
                end
                updateplot=updateplot&(~isempty(obj.dataObject) || (isfield(obj.dataObject,'status') && obj.dataObject.status));
                obj.draw3Dmovie=updateplot;
                
                if ~updateplot
                    obj.deletetry(obj.handels.selectedSpotTab);
                end
                
                if ischar(answer{1})
                    ws = load(answer{1});
                    ws=ws.vs;
                elseif isstruct(src)
                    ws=src;
                end
    
                if ~strcmpi(ws.class,obj.classSingle) && ~strcmpi(ws.class, obj.classAdded) 
                    error(['This is no ' class ' class'])
                else
%                     obj.class = ws.class;
                    obj.startFrame = ws.startFrame;
                    obj.PSFSigma = ws.PSFSigma;
                    obj.detectionFrame = ws.detectionFrame;
                    obj.NFrames = ws.NFrames;

                    obj.AddedTrackingInfo.spots(length(obj.AddedTrackingInfo.files)) = length(ws.spotsIncluded);
                    
                    obj.flipCameras = ws.flipCameras;
                    obj.subRegionSelection = ws.subRegionSelection;
                    obj.colocalized = ws.colocalized;
                    obj.rastergramStartframe=ws.rastergramStartframe;
                    
                    obj.numberOfSpotsForRastergram=ws.numberOfSpotsForRastergram;
                    obj.THP=ws.THP;
                    obj.THX=ws.THX;
                    obj.THBG=ws.THBG;
                    obj.THS=ws.THS;
                    obj.Nbootstraps=ws.Nbootstraps;
                    obj.tagetTrace =ws.tagetTrace;
                    obj.colorTarget = ws.colorTarget;
                    obj.colorComplex= ws.colorComplex;


                    %results
                    obj.sortMap=ws.sortMap;
                    obj.spotNrs=ws.spotNrs;
                    obj.meanTime2FirstEvent=ws.meanTime2FirstEvent;
                    obj.meanDwellTime=ws.meanDwellTime;
                    obj.rawFitResultsCam1=ws.rawFitResultsCam1;
                    obj.rawFitResultsCam2=ws.rawFitResultsCam2;
                    obj.coordsCam1=ws.coordsCam1;
                    obj.coordsCam2=ws.coordsCam2;
                    obj.THPmax=ws.THPmax;
                    obj.THXmax=ws.THXmax;
                    obj.THSmax=ws.THSmax;
                    obj.THBGmax=ws.THBGmax;


                    obj.coord1=ws.coord1;
                    obj.coord2=ws.coord2;
                    obj.spotsIncluded=ws.spotsIncluded;
                    obj.maskFilt2=ws.maskFilt2;
                    obj.maskFilt1=ws.maskFilt1;
                    obj.photons1I=ws.photons1I;
                    obj.photons2I=ws.photons2I;
                    obj.bg2I=ws.bg2I;
                    obj.bg1I=ws.bg1I;
                    obj.Sigma1I=ws.Sigma1I;
                    obj.deltax=ws.deltax;
                    obj.deltay=ws.deltay;
                    obj.seq=ws.seq;
                    obj.intergrationTime=ws.intergrationTime;
                    obj.segmentedSpots=ws.segmentedSpots;
                    obj.coordsComplex=ws.coordsComplex;
                    obj.coordsTarget=ws.coordsTarget;
                    obj.handels.selectedSpotIndex = [];


                    maxNumberOfImages=size(obj.spotsIncluded,1);
                    set(obj.handels.spotSelectionSlider, 'SliderStep', [1/maxNumberOfImages 10/maxNumberOfImages]);       
                    set(obj.handels.spotSelectionSlider, 'Max', maxNumberOfImages);       
                    set(obj.handels.spotSelectionSlider, 'Value', round(maxNumberOfImages/2)); 
%                         obj.imageStackCam2=ws.imageStackCam2;
%                         obj.imageStackCam1=ws.imageStackCam1;

                    obj.handels.Figure.UserData.sliderEventClass.changeLimits([obj.THPmax; obj.THXmax; obj.THBGmax; obj.THSmax]);
                    obj.handels.Figure.UserData.sliderEventClass.changePosition([obj.THP/(obj.THPmax-1)*get(obj.handels.intensityTHSlider,'Max');...
                    obj.THX/(obj.THXmax-1)*get(obj.handels.positionTHSlider,'Max');...
                    obj.THBG/(obj.THBGmax-1)*get(obj.handels.backgroundRatioTHSlider,'Max');...
                    obj.THS/(obj.THSmax-1)*get(obj.handels.sigmaTHSlider,'Max')]);
    %                 if updateplot
                    obj.plotImages;
                    obj.status=1;
    %                 end
                end
            else
                ws=[];
            end
        end
        
        function addvars(obj,src,updateplot)
            if ~obj.addedvars
                class = obj.classSingle;
            else
                class = obj.classAdded;
            end
                if ishandle(src)

                [file, pathname] = uigetfile([class '.mat'] );
                if file ~= 0
                    answer{1} = fullfile(pathname, file);
                else
                    answer=[];
                end
                
                elseif ischar(src)
                    answer{1} = src;
                else
                    answer{1} = [];  
                end

                if ~isempty(answer)
                    if isfield(obj.AddedTrackingInfo,'files') && size(obj.AddedTrackingInfo.files,1) > 0
                        obj.AddedTrackingInfo.files{end+1} = answer;
                    else
                        obj.AddedTrackingInfo.files = [];
                        obj.AddedTrackingInfo.files{1} = answer;
                    end
                    updateplot=false;
                    if nargin < 3 || isempty(updateplot) || isobject(updateplot)
                         obj.draw3Dmovie=true;
                     else
                         obj.draw3Dmovie=updateplot;
                    end

                    if ischar(answer{1})
                        ws = load(answer{1});
                        ws=ws.vs;
                    elseif isstruct(src)
                        ws=src;
                    end
                    if ~strcmpi(ws.class,obj.classSingle) && ~strcmpi(ws.class, obj.classAdded) 
                        error(['This is no ' class ' class'])
                    else
                        
                        obj.AddedTrackingInfo.spots(length(obj.AddedTrackingInfo.files)) = length(ws.spotsIncluded);
                        obj.rawFitResultsCam1=appendRawFitResults(obj.rawFitResultsCam1,ws.rawFitResultsCam1);
                        obj.rawFitResultsCam2=appendRawFitResults(obj.rawFitResultsCam2,ws.rawFitResultsCam2);
                        obj.coordsCam1=cat(1,obj.coordsCam1,ws.coordsCam1);
                        obj.coordsCam2=cat(1,obj.coordsCam2,ws.coordsCam2);

                        obj.coord1=cat(1,obj.coord1,ws.coord1);
                        obj.coord2=cat(1,obj.coord2,ws.coord2);
                        obj.spotsIncluded=cat(1,obj.spotsIncluded,ws.spotsIncluded);
                        obj.maskFilt2=cat(1,obj.maskFilt2,ws.maskFilt2);
                        obj.maskFilt1=cat(1,obj.maskFilt1,ws.maskFilt1);
                        obj.photons1I=cat(1,obj.photons1I,ws.photons1I);
                        obj.photons2I=cat(1,obj.photons2I,ws.photons2I);
                        obj.bg2I=cat(1,obj.bg2I,ws.bg2I);
                        obj.bg1I=cat(1,obj.bg1I,ws.bg1I);
                        obj.Sigma1I=cat(1,obj.Sigma1I,ws.Sigma1I);
                        obj.deltax=cat(1,obj.deltax,ws.deltax);
                        obj.deltay=cat(1,obj.deltay,ws.deltay);
                        obj.seq=cat(1,obj.seq,ws.seq);
                        obj.rastergramStartframe = cat(2,obj.rastergramStartframe,ws.rastergramStartframe);

                        obj.coordsComplex=cat(1,obj.coordsComplex,ws.coordsComplex);
                        obj.coordsTarget=cat(1,obj.coordsTarget,ws.coordsTarget);
                        obj.handels.selectedSpotIndex = [];

                        obj.plotImages;
                        obj.status=1;
                        obj.addedvars = true;
                        
                        obj.deletetry(obj.handels.selectedSpotTab);
                        obj.settry(obj.handels.childMenu(5),'Separator','on');
                        obj.deletetry(obj.handels.childMenu(1));
                        obj.deletetry(obj.handels.childMenu(4));
                    end 
                end
                
            function res = appendRawFitResults(rawFitResultsCamBase,rawFitResultsCamNew)
                res.RoiStart = cat(1,rawFitResultsCamBase.RoiStart,rawFitResultsCamNew.RoiStart);
                res.ROIStack = cat(3,rawFitResultsCamBase.ROIStack,rawFitResultsCamNew.ROIStack);
                res.Coord = cat(1,rawFitResultsCamBase.Coord,rawFitResultsCamNew.Coord);
                res.Photons = cat(1,rawFitResultsCamBase.Photons,rawFitResultsCamNew.Photons);
                res.Bg = cat(1,rawFitResultsCamBase.Bg,rawFitResultsCamNew.Bg);
                res.Sigma = cat(1,rawFitResultsCamBase.Sigma,rawFitResultsCamNew.Sigma);
                res.Sigma_STD = cat(1,rawFitResultsCamBase.Sigma_STD,rawFitResultsCamNew.Sigma_STD);
                res.Frame = cat(1,rawFitResultsCamBase.Frame,rawFitResultsCamNew.Frame);
                res.CRLB_STD = cat(1,rawFitResultsCamBase.CRLB_STD,rawFitResultsCamNew.CRLB_STD);
                res.Photons_STD = cat(1,rawFitResultsCamBase.Photons_STD,rawFitResultsCamNew.Photons_STD);
                res.Bg_STD = cat(1,rawFitResultsCamBase.Bg_STD,rawFitResultsCamNew.Bg_STD);
                res.LL = cat(2,rawFitResultsCamBase.LL,rawFitResultsCamNew.LL);

            end
        end
        function exportData(obj,src,~)
            traces2_allfr = obj.getTraceMatrix;
            
            filename = uiputfile('*.mat')
            if ~isempty(filename)
                data.Trace = traces2_allfr;
                data.Background = obj.bg1I(obj.spotsIncluded,:,2-0)*(obj.PSFSigma*2)^2;
                data.Intensity = obj.photons1I(obj.spotsIncluded,:,2-0);
                data.DxSq =  obj.deltay(obj.spotsIncluded,:).^2+obj.deltax(obj.spotsIncluded,:).^2;
                data.Sigma =  obj.Sigma1I(obj.spotsIncluded,:,2-0);
                data.spotIDs = obj.spotNrs;
                save(filename,'data')
            end
        end

        function exportRastergramFig(obj,src,~)
            filename = uiputfile('*.pdf')
            if ~isequal(filename,0)
                h3 = figure;
                h1=obj.handels.rastergramPlot.Parent;
%                 objects=allchild(h1);
                copyobj(get(h1,'children'),h3);
                idx=[];
                for i = 1:length(h3.Children)
                    if isa(h3.Children(i),'matlab.ui.control.UIControl')
                       idx=[idx,i];
                    end
                end
                delete(h3.Children(idx))
                
                export_fig(filename,'-pdf','-transparent', h3)
%                 export_fig(filename,'-eps','-transparent', obj.handels.rastergramPlot)
            end
        end
        
        function exportRastergramData(obj,src,~)
            filename = uiputfile('*.mat')
            if ~isequal(filename,0)
                data.plotData = obj.handels.rastergramPlot.UserData.plotData;
                save(filename,'data')
            end
        end
        
        function exportTraceFigure(obj,src,~)
            filename = uiputfile('*.pdf')
            if ~isequal(filename,0)
                h3 = figure;
                h1=obj.handels.traceEstimateTabPlot.Parent;
%                 objects=allchild(h1);
                copyobj(get(h1,'children'),h3);
                idx=[];
                for i = 1:length(h3.Children)
                    if isa(h3.Children(i),'matlab.ui.control.UIControl')
                       idx=[idx,i];
                    end
                end
                delete(h3.Children(idx))
                
                export_fig(filename,'-pdf','-transparent', h3)
            end
        end
        
        function exportTraceData(obj,src,~)
            filename = uiputfile('*.mat')
            if ~isequal(filename,0)
                traces2_allfr = obj.getBinary(obj.handels.selectedSpotIndex);
                data.Trace = traces2_allfr;
                data.Background = obj.bg1I(obj.spotsIncluded,:,2-0)*(obj.PSFSigma*2)^2;
                data.Intensity = obj.photons1I(obj.spotsIncluded,:,2-0);
                data.DxSq =  obj.deltay(obj.spotsIncluded,:).^2+obj.deltax(obj.spotsIncluded,:).^2;
                data.Sigma =  obj.Sigma1I(obj.spotsIncluded,:,2-0);
                data.spotIDs = obj.spotNrs;
                save(filename,'data')
            end
        end
            
    end
end
