classdef guiDarkROIs < handle
   properties
    Figure
    Axes
    CLimit
    dataObject
    htab1
    selRec
    startFrame
    detectionFrame
    hp
    colocalized
    sortMap
    spotNrs
    hras
    meanTime2FirstEvent
    meanDwellTime
    spotLine
    
    NFrames
    himg1
    himg2
    C
    alignmentObject
    driftObject
    spotDetectorObject
    rawFitResultsCam1
    rawFitResultsCam2
    coordsCam1
    coordsCam2
    lh
    lh2
    Nspot
    hh
    hAxis2
    htab2
    coord1
    hsl1
    s
    s1
    s2
    THPmax
    THXmax
    THSmax
    THBGmax
    flip
    
    labs
    cbh
    htab3
    spotsIncluded
    timeLine
    initialClick
    amin1
    amax1
    
    amin2
    amax2    
    %tajectory variables 
    maskFilt2
    maskFilt1
    intergrationTime
    h1
    h2
    
    photons1I
    photons2I

    bg2I
    bg1I

    deltax
    deltay

    seq
   end
   
	events
       ParamChange
       timePointChanged
       spotChanged
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
       
       function selectSubRegion(obj,C,~)
            if size(C,1) == 2 && size(C,2) == 2
                obj.C = C;
            else
                h=dipshow(joinchannels('RGB',squeeze(obj.alignmentObject.rgbImage(:,:,1)),squeeze(obj.alignmentObject.rgbImage(:,:,2))),'lin');
                [~,obj.C] = dipcrop(h);
                close(h)
            end

       end
       
        function  setParamsGeneral(obj,src,~)
            prompt={'Enter a value of max I','Enter a value of max \sigma','Enter a value of max \Delta X'};
            name = 'Threshold Params';
                

            defaultans = {num2str(obj.THPmax), num2str(obj.THXmax), num2str(obj.THSmax)}; %1e3,10, 5
            options.Interpreter = 'tex';
            
            answer = inputdlg(prompt,name,[1 40],defaultans,options);
            obj.THPmax = str2num(answer{1});
            obj.THXmax = str2num(answer{2});
            obj.THSmax = str2num(answer{3});
        end

       
        function obj = guiDarkROIs(h,dataObject,alignmentObject,driftObject,spotDetectorObject,intergrationTime,menuhandle,menuName,startFrame,detectionFrame,colocalized,flip)
            
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
            obj.colocalized=colocalized;
            obj.flip=flip;
            obj.intergrationTime=intergrationTime;
            obj.startFrame=startFrame;
            obj.detectionFrame=detectionFrame;
            
            obj.alignmentObject=alignmentObject;
            obj.driftObject=driftObject;
            if isa(spotDetectorObject,'spotDetectorMultiplexRNAWrapper')
                spotDetectorObject = spotDetectorObject.spotDetectorObjects{end};
            end
            obj.spotDetectorObject=spotDetectorObject;
            
            
            obj.THPmax=1e3;
            obj.THXmax=10;
            obj.THSmax=5;
            obj.THBGmax=5;
            
            obj.dataObject = dataObject;
            topTab = uitab(h, 'Title', ['Analyze Dark ' menuName]);
            htabgroup = uitabgroup(topTab);
            obj.htab1 = uitab(htabgroup, 'Title', 'Estimated Trace');
            obj.htab2 = uitab(htabgroup, 'Title', 'Selected Spot');
            obj.htab3 = uitab(htabgroup, 'Title', 'Rastergram');
            
            if nargin < 7 || isempty(menuhandle)
                f = uimenu('Label','Analyze Dark ');
            else
                f = uimenu('Parent',menuhandle,'Label',menuName);
            end
            uimenu('Parent',f,'Label','Select Subregion','Callback',@obj.selectSubRegion);
            uimenu('Parent',f,'Label','Threshold params','Callback',@obj.setParamsGeneral);
            uimenu('Parent',f,'Label','Analyze','Callback',@obj.analyze);
            uimenu('Parent',f,'Label','Calculate Rastergram','Callback',@obj.calcras);

            
            obj.lh = addlistener(obj,'timePointChanged',@obj.plotTimePoint)
            obj.lh2 = addlistener(obj,'spotChanged',@obj.plotImages)
            
            slmin = 1;
            slmax = 100;
            obj.hsl1 = uicontrol('Parent', obj.htab1,'Style','slider','Min',slmin,'Max',slmax,...
                            'SliderStep',[1 1]./(slmax-slmin),'Value',51,...
                            'Units','normalized','Position',[ 0.1 0.01 0.8 0.05]);
     
            slmin = 0;
            slmax = 100;
            obj.s = uicontrol('Parent',  obj.htab1,'Style','Slider','SliderStep',[1 1]./(slmax-slmin),...
            'Units','normalized','Position',[0 0 0.01 0.8],'Min',slmin,'Max',slmax,...
                  'Value',51,'String','a');                           
            obj.s1 = uicontrol('Parent',  obj.htab1,'Style','Slider','SliderStep',[1 1]./(slmax-slmin),...
            'Units','normalized','Position',[0.01 0 0.01 0.8],'Min',slmin,'Max',slmax,...
                  'Value',51); 
            obj.s2 = uicontrol('Parent',  obj.htab1,'Style','Slider','SliderStep',[1 1]./(slmax-slmin),...
            'Units','normalized','Position',[0.02 0 0.01 0.8],'Min',slmin,'Max',slmax,...
                  'Value',0); 

            THP = (obj.THPmax-1)/get(obj.s,'Max')*get(obj.s,'Value');
            THX = (obj.THXmax-1)/get(obj.s,'Max')*get(obj.s1,'Value');              
            obj.labs = uicontrol('Parent',  obj.htab1,'Style','text','Units','normalized','position',[ 0 0.8 0.05 0.1],'String',sprintf('I=%0.2g\n X=%0.2g',THP,THX));
            obj.cbh = uicontrol('Parent',  obj.htab1,'Style','checkbox',... 
                'String','Include trace','Units','normalized',... 
                'Value',0,'Position',[0 .9 0.1 0.1]);
            
            set(obj.cbh,'Callback',@obj.includeSpot)
            
            set(obj.hsl1,'Callback',@obj.changeSpot)
            set(obj.s,'Callback',@obj.plotImages)
            set(obj.s1,'Callback',@obj.plotImages)
            set(obj.s2,'Callback',@obj.plotImages)
            
            obj.Figure = getParentFigure(h);
        end
        
        function includeSpot(obj,src,~)
            if ~isempty(obj.spotsIncluded) && ~isempty(obj.Nspot)
                obj.spotsIncluded(obj.Nspot) = get(obj.cbh,'Value');
                if obj.spotsIncluded(obj.Nspot) == 1
                    line(obj.coord1(obj.Nspot,1),obj.coord1(obj.Nspot,2),...
                    'Color','y','marker','o','LineStyle','none','Parent',obj.hAxis2,'markersize',4,'tag','pointsSM','buttondownfcn',@obj.selectMarker)
                elseif obj.spotsIncluded(obj.Nspot) == 0
                    line(obj.coord1(obj.Nspot,1),obj.coord1(obj.Nspot,2),...
                    'Color','r','marker','o','LineStyle','none','Parent',obj.hAxis2,'markersize',4,'tag','pointsSM','buttondownfcn',@obj.selectMarker)
                end
            end
        end
        
        function changeSpot(obj,src,~)
            if  ~isempty(squeeze(obj.coord1))
                 pos = round(get(src,'Value'));
                obj.Nspot = round((size(obj.coord1,1))/get(src,'Max')*pos);
                notify(obj,'spotChanged');
            end
        end
        
        function dataChange(obj,src,~)
            notify(obj,'ParamChange');
        end

        
         function  savevars(obj,src,~)
        end
       

      
        function  analyze(obj,src,~)
            obj.clear;
            coords =   obj.spotDetectorObject.darkCoords;
            frames =   (obj.startFrame-1).*ones(size(coords,1),1);
            coords_1st=coords(frames==obj.startFrame-1,:);
     

            tSegmentedSpots = transformPointsInverse(obj.alignmentObject.tformTotal,coords_1st);
            if ~isempty(obj.C)
                idx = obj.C(1,1) < tSegmentedSpots(:,1) &...
                obj.C(1,1)+obj.C(2,1) > tSegmentedSpots(:,1)&...
                obj.C(1,2) < tSegmentedSpots(:,2) &...
                obj.C(1,2)+obj.C(2,2) > tSegmentedSpots(:,2);
            else
                idx = true(size(tSegmentedSpots,1),1);
            end
            maxNumberOfImages=sum(idx);
            obj.spotsIncluded = true(maxNumberOfImages,1);
            set(obj.hsl1, 'SliderStep', [1/maxNumberOfImages 10/maxNumberOfImages]);       
            set(obj.hsl1, 'Max', maxNumberOfImages);       
            set(obj.hsl1, 'Value', round(maxNumberOfImages/2)); 
            
            segmentedSpots = coords_1st(idx,:);

            %%
            obj.coordsCam1=[];
            obj.coordsCam2=[];
            for w=obj.startFrame:size(obj.dataObject.dataCam{1},3)
                coords_nth = segmentedSpots-repmat(obj.driftObject.drift(w,:),size(segmentedSpots,1),1);
                obj.coordsCam1=cat(1,obj.coordsCam1,[coords_nth (w-1).*ones(size(segmentedSpots,1),1)]);
                obj.coordsCam2=cat(1,obj.coordsCam2,[transformPointsInverse(obj.alignmentObject.tformTotal,coords_nth) (w-1).*ones(size(segmentedSpots,1),1)]); 
            end

            % MLE Fit Intensities Cam1
            paramsFit = obj.spotDetectorObject.paramsFit;
            paramsFit.FitSigma=true;
             coords1 = [obj.coordsCam1(:,1) obj.coordsCam1(:,2) obj.coordsCam1(:,3)];
            [ obj.rawFitResultsCam1 ] = fitBoxCenters( single(squeeze(obj.dataObject.dataCam{1})),coords1,paramsFit);
            
             coords2 = [obj.coordsCam2(:,1) obj.coordsCam2(:,2) obj.coordsCam2(:,3)];
            [ obj.rawFitResultsCam2 ] = fitBoxCenters( single(squeeze(obj.dataObject.dataCam{2})),coords2,paramsFit);
        
            obj.NFrames  =obj.rawFitResultsCam2.Frame(end)+2-obj.startFrame;
            obj.maskFilt2 = true(size(obj.rawFitResultsCam2.Photons,1),1);
            obj.maskFilt1 = true(size(obj.rawFitResultsCam1.Photons,1),1);

            photons1 = obj.rawFitResultsCam1.Photons(obj.maskFilt1,:);
            bg1 = obj.rawFitResultsCam1.Bg(obj.maskFilt1,:);
            bg2 = obj.rawFitResultsCam2.Bg(obj.maskFilt2,:);
            photons2 = obj.rawFitResultsCam2.Photons(obj.maskFilt2,:);

            obj.photons1I(:,:,1+obj.flip) = reshape(photons1,[size(photons1,1)/obj.NFrames obj.NFrames]);
            obj.photons1I(:,:,2-obj.flip) = reshape(photons2,[size(photons1,1)/obj.NFrames obj.NFrames]);

            Sigma = reshape(sqrt(sum(obj.rawFitResultsCam2.Sigma.^2,2)),[size(photons1,1)/obj.NFrames obj.NFrames]);
            obj.bg1I(:,:,2-obj.flip) = reshape(bg2,[size(photons2,1)/obj.NFrames obj.NFrames]);
            obj.bg1I(:,:,1+obj.flip) = reshape(bg1,[size(photons2,1)/obj.NFrames obj.NFrames]);
            
            
            if obj.colocalized
                delta = transformPointsInverse(obj.alignmentObject.tformTotal,obj.rawFitResultsCam1.Coord)-obj.rawFitResultsCam2.Coord;
            elseif ~obj.flip && ~obj.colocalized 
                delta = obj.rawFitResultsCam1.Coord-coords1(:,1:2);
            elseif obj.flip && ~obj.colocalized
                 delta = obj.rawFitResultsCam2.Coord-coords2(:,1:2);
            end
            
            obj.deltax = reshape(delta(:,1),[size(delta(:,1),1)/obj.NFrames obj.NFrames]);
            obj.deltay = reshape(delta(:,2),[size(delta(:,2),1)/obj.NFrames obj.NFrames]);

            obj.seq = createSeq(obj.rawFitResultsCam1,obj.rawFitResultsCam2,obj.coordsCam2,obj.alignmentObject.tformTotal,obj.startFrame-1);
                      
            
            obj.plotImages;

        end
        
        function rs = getSummary(obj)
        end
        
        function clear(obj,src,~)
            obj.rawFitResultsCam1 = [];
            obj.rawFitResultsCam2 = [];
            obj.coordsCam1 = [];
            obj.coordsCam2 = [];
            obj.Nspot = [];
            obj.coord1  = [];
            obj.spotsIncluded = [];

        end
        
        function plotImages(obj,src,~)
            reDraw = isempty(obj.Nspot);

            if reDraw
                obj.Nspot =  get(obj.hsl1,'Value');
                maxSpots =  get(obj.hsl1,'Max');
                spotIdx = obj.Nspot:maxSpots:obj.NFrames*maxSpots;
                obj.amin1 = min(dip_image(obj.rawFitResultsCam1.ROIStack(:,:,spotIdx)));
                obj.amax1 = max(dip_image(obj.rawFitResultsCam1.ROIStack(:,:,spotIdx)));

                obj.amin2 = min(dip_image(obj.rawFitResultsCam2.ROIStack(:,:,spotIdx)));
                obj.amax2 = max(dip_image(obj.rawFitResultsCam2.ROIStack(:,:,spotIdx)));
            end
            
            set(obj.cbh,'Value',obj.spotsIncluded(obj.Nspot))
            THP = (obj.THPmax-1)/get(obj.s,'Max')*get(obj.s,'Value');
            THX = (obj.THXmax-1)/get(obj.s1,'Max')*get(obj.s1,'Value');
            THBG = (obj.THBGmax-1)/get(obj.s2,'Max')*get(obj.s2,'Value');
            set(obj.labs,'String',sprintf('I=%0.2g, X=%0.2g, Ph x BG, %0.2g',THP,THX,THBG))
            
            traces2_allfr = obj.photons1I(obj.Nspot,:,2-0) > THP...
                & sqrt(obj.deltay(obj.Nspot,:).^2+obj.deltax(obj.Nspot,:).^2) < THX...
                & obj.photons1I(obj.Nspot,:,2-0) > THBG*mean(obj.bg1I(obj.Nspot,:,2-0))*(obj.spotDetectorObject.paramsFit.PSFSigma*2)^2; %...

            ymax = max(max(obj.bg1I(obj.Nspot,:,2-0)*(obj.spotDetectorObject.paramsFit.PSFSigma*2)^2),max(max(obj.photons1I(obj.Nspot,:,2-0)),max(obj.photons1I(obj.Nspot,:,1+0))));
            delete(obj.hh);
            obj.hh = subplot(2,2,1:2,'Parent', obj.htab1,'buttondownfcn',@obj.clickline,'xlimmode','manual');
            
            axis([0 size(obj.photons1I,2) 0 ymax])

            line(1:size(obj.photons1I,2),obj.photons1I(obj.Nspot,:,2-0),'Color','r','Marker','x','parent',obj.hh)
            hold on
            line(1:size(obj.photons1I,2),obj.photons1I(obj.Nspot,:,1+0),'Color','g','parent',obj.hh,'tag','background','hittest','off')
            line(1:size(obj.photons1I,2),obj.bg1I(obj.Nspot,:,2-0)*(obj.spotDetectorObject.paramsFit.PSFSigma*2)^2,'Color',[.7 .5 0],'parent',obj.hh,'tag','background','hittest','off')
            line(1:size(obj.photons1I,2),traces2_allfr*ymax,'Color','b','parent',obj.hh,'tag','background','hittest','off')

            if ~reDraw & ishandle(obj.timeLine)
               xcoords = get(obj.timeLine,'xdata');
               ii= xcoords(1);
            else            
                ii=obj.startFrame+10;
            end
            
            obj.timeLine = line([ii ii],[0 ymax],'LineStyle','--','Color','k','Linewidth',1,'parent',obj.hh,'tag','timeLine','hittest','off');
            ylabel('I [# Photons]')
            xlabel('Time [s]')
            drawnow
  
             legend('Component','Reference','Background','Binary','Orientation','horizontal') 

            if size(obj.seq,3) > 3
                b = sqrt(obj.seq(obj.Nspot,:,4).^2+obj.seq(obj.Nspot,:,5).^2);
            end
            hold off
         
            if isempty(obj.hAxis2) || reDraw
                idx = obj.coordsCam1(:,3) == obj.startFrame;
                obj.hAxis2 = subplot(1,1,1,'Parent', obj.htab2,'buttondownfcn',@obj.selectMarker,'xlimmode','manual')
                hh = subimage(imadjust(mat2gray(double(stretch(obj.dataObject.dataCam{1}(:,:,obj.detectionFrame))))))  
                set(hh,'buttondownfcn',@obj.selectMarker) 
                obj.coord1(:,1) = obj.coordsCam1(idx,1);
                obj.coord1(:,2)=obj.coordsCam1(idx,2);

                line(obj.coord1(:,1),obj.coord1(:,2),...
                    'Color','y','marker','o','Parent',obj.hAxis2,'LineStyle','none','markersize',4,'tag','pointsSM','buttondownfcn',@obj.selectMarker)
                x=obj.coord1(obj.Nspot,1)
                y=obj.coord1(obj.Nspot,2)
                boxsize=10;
                obj.selRec = rectangle('Position',[x-boxsize/2 y-boxsize/2 boxsize boxsize],'EdgeColor','y','tag','selectionBox','hittest','off','buttondownfcn',@obj.selectMarker)
            else
                boxsize=10;
                set(obj.selRec,'Position',[obj.coord1(obj.Nspot,1)-boxsize/2 obj.coord1(obj.Nspot,2)-boxsize/2 boxsize boxsize]);
            end
            obj.plotTimePoint;  
        end      

        function selectMarker(obj,src,ev)
            clicked=get(gca,'currentpoint');
            xcoord=clicked(1,1,1);
            ycoord=clicked(1,2,1);

            obj.Nspot = knnsearch(obj.coord1,[xcoord ycoord]) 
            boxsize=10;
            set(obj.selRec,'Position',[obj.coord1(obj.Nspot,1)-boxsize/2 obj.coord1(obj.Nspot,2)-boxsize/2 boxsize boxsize]);

            notify(obj,'spotChanged');
          end
        
        function clickline(obj,src,ev)
            clicked=get(gca,'currentpoint');
            obj.initialClick=clicked(1,1,1);
           
            set(gcf,'windowbuttonmotionfcn',@obj.dragline)
            set(gcf,'windowbuttonupfcn',@obj.dragdone)
        end
        
        function dragline(obj,src,ev)          
            clicked=get(gca,'currentpoint');
            xcoord=clicked(1,1,1);
            
            set(obj.timeLine,'xdata',[xcoord xcoord]);
            if size(obj.seq,3) > 3
                ii = round(xcoord);
                b = sqrt(obj.seq(obj.Nspot,ii,4).^2+obj.seq(obj.Nspot,ii,5).^2);
                title(sprintf('spot ID = %d\n Frame = %d, \\Delta X = %0.2g, I = %0.2g, bg=%0.2g',obj.Nspot,ii,b,obj.photons1I(obj.Nspot,ii,2-0),obj.bg1I(obj.Nspot,ii,2-0)*(obj.spotDetectorObject.paramsFit.PSFSigma*2)^2),'Parent',obj.hh);
            end            
        end


        function dragdone(obj,src,ev)
            clicked=get(gca,'currentpoint');
            finalClick=clicked(1,1,1);
            
            if finalClick == obj.initialClick
                xcoord = get(obj.timeLine,'xdata')
                if finalClick < xcoord(1)
                    set(obj.timeLine,'xdata',xcoord-1);
                else
                    set(obj.timeLine,'xdata',xcoord+1);
                end
            end

            set(gcf,'windowbuttonmotionfcn','');
            set(gcf,'windowbuttonupfcn','');
            notify(obj,'timePointChanged');
        end
        
        function calcras(obj,src,~)
            THP = (obj.THPmax-1)/get(obj.s,'Max')*get(obj.s,'Value');
            THX = (obj.THXmax-1)/get(obj.s1,'Max')*get(obj.s1,'Value');
            THBG = (obj.THBGmax-1)/get(obj.s2,'Max')*get(obj.s2,'Value');
            
             traces2_allfr = obj.photons1I(obj.spotsIncluded,:,2-0) > THP...
                & sqrt(obj.deltay(obj.spotsIncluded,:).^2+obj.deltax(obj.spotsIncluded,:).^2) < THX...
                & obj.photons1I(obj.spotsIncluded,:,2-0) > repmat(THBG*mean(obj.bg1I(obj.spotsIncluded,:,2-0),2)*(obj.spotDetectorObject.paramsFit.PSFSigma*2)^2,[1 size(obj.bg1I,2) 1]); %...

            obj.spotNrs = find(obj.spotsIncluded);
            
            if isempty(obj.intergrationTime)
                obj.intergrationTime=1;
            elseif ischar(obj.intergrationTime)
                ws = load(obj.intergrationTime,'vid')
                ttb = ws.vid.ttb;
            elseif isnumeric(obj.intergrationTime)
                ttb=1000*[0:obj.intergrationTime:obj.intergrationTime*(size(traces2_allfr,2)-1)]; % This is an array containing relative time stamps for all frames
            end
                           
            nframes=length(ttb); 

            % Code provided by the Jeff Gelles Lab:
            traces2_allfr(:,end)=[]; % because time refers to the start of the frame, need to do this so that each frame has start time and end time 
            traces2_allfr_s=num2str(traces2_allfr); 
            traces2_allfr_sc=num2cell(traces2_allfr_s,2); % cell array is needed for regexp to work
            traces2_allfr_sc=regexprep(traces2_allfr_sc,' ',''); % this makes a searchable cell array where each cell is a binary trace in string format

            Intervals.CumulativeIntervalArray=[]; %initialize array

            for i=1:length(traces2_allfr_sc)
            % Key (Larry's flags):
            % "1..." = -3
            % "...1" = 3
            % "0..." = -2
            % "...0" = 2
            % "1...0...1" = 0
            % "0...1...0" = 1
            exp=traces2_allfr_sc{i,1};
            [start_flag_3, end_flag_3]=regexp(exp,'^1+');
            [start_flag_2, end_flag_2]=regexp(exp,'^0+');
            [start_flag0, end_flag0]=regexp(exp,'(?<=1)0{2,}(?=1)'); % at least 2 zeros flanked by ones, 1 frame events are discarded
            [start_flag1, end_flag1]=regexp(exp,'(?<=0)1{2,}(?=0)'); % at least 2 ones flanked by zeros, 1 frame events are discarded
            [start_flag2, end_flag2]=regexp(exp,'0+$');
            [start_flag3, end_flag3]=regexp(exp,'1+$');
            if ~isempty(start_flag_3)
                Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [-3 start_flag_3 end_flag_3 ...
                    end_flag_3-start_flag_3+1 ttb(end_flag_3+1)-ttb(start_flag_3) 0 i]];
            end
            if ~isempty(start_flag_2)
                Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [-2 start_flag_2 end_flag_2 ...
                    end_flag_2-start_flag_2+1 ttb(end_flag_2+1)-ttb(start_flag_2) 0 i]];
            end
            if ~isempty(start_flag0)
                Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [zeros(length(start_flag0),1) ...
                    start_flag0' end_flag0' end_flag0'-start_flag0'+1 ttb(end_flag0+1)'-ttb(start_flag0)' ...
                    zeros(length(start_flag0),1) i*ones(length(start_flag0),1)]];
            end
            if ~isempty(start_flag1)
                Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [ones(length(start_flag1),1) ...
                    start_flag1' end_flag1' end_flag1'-start_flag1'+1 ttb(end_flag1+1)'-ttb(start_flag1)' ...
                    zeros(length(start_flag1),1) i*ones(length(start_flag1),1)]];
            end
            if ~isempty(start_flag2) & start_flag2~=start_flag_2
                Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [2 start_flag2 end_flag2 ...
                    end_flag2-start_flag2+1 ttb(end_flag2+1)-ttb(start_flag2) 0 i]];
            end
            if ~isempty(start_flag3) & start_flag3~=start_flag_3
                Intervals.CumulativeIntervalArray=[Intervals.CumulativeIntervalArray; [3 start_flag3 end_flag3 ...
                    end_flag3-start_flag3+1 ttb(end_flag3+1)-ttb(start_flag3) 0 i]];
            end
            end

            Intervals.CumulativeIntervalArray=sortrows(Intervals.CumulativeIntervalArray,[7,2]);

            cia=Intervals.CumulativeIntervalArray;
            logik=[];
            logik=(cia(:,1)==-2);
            Time2FirstEvent=cia(logik,5); 
            time_empty=max(Time2FirstEvent)*ones(length(Time2FirstEvent),1);
            t=(Time2FirstEvent~=time_empty);
            time21st=Time2FirstEvent(t);
            obj.meanTime2FirstEvent = mean(time21st);


            cia=Intervals.CumulativeIntervalArray;
            logik=[];
            logik=cia(:,1)==1; % &cia(:,4)>1;
            DwellTime=cia(logik,5);
            obj.meanDwellTime = mean(DwellTime);


            logik1st=cia(:,1)==-2; % this is to sort in order of first arrival time
            Arrivals=[cia(logik1st,3) cia(logik1st,7)]; %1st column is arrival frame, 2nd column is aoi#
            [Sorted_1stevent,obj.sortMap]=sortrows(Arrivals,1); % Sort in the order of  arrival time for 1st arrival
            newcia=[]; % initialize new cumulative interval array for sorted IDS
            for i=1:max(size(Sorted_1stevent));
                logik=cia(:,7)==Sorted_1stevent(i,2);
                instant_cia=cia(logik,:); instant_cia(:,7)=i;
                newcia=[newcia;instant_cia]; 
            end
            Intervals.CumulativeIntervalArray=newcia;

            vid.ttb=ttb; %1000*[0:1.5:898.5]; % This is an array containing relative time stamps for all frames
            vid.nframes=length(vid.ttb); % This is total number of frames
            save swaprandvidfile.mat vid

            num_traces=max(Intervals.CumulativeIntervalArray(:,7));
            save IDS2_sorted.mat Intervals;
            obj.hras = subplot(1,1,1,'Parent', obj.htab3,'buttondownfcn',@obj.clickhline,'xlimmode','manual');
            func_graph_intervals_nooverlap_v1pt0(obj.hras,'IDS2_sorted','swaprandvidfile.mat',1,vid.nframes,num_traces,[1 1 1],[1 1 1],[0 0 0],0); % make rastergram
            obj.spotLine = line([0 size(obj.photons1I,2)],[1 1].*find(obj.Nspot == obj.sortMap),'LineStyle','--','Color','r','Linewidth',1,'parent',obj.hras,'tag','timeLine','hittest','off');
            title(sprintf('spot ID = %d, meanTime2FirstEvent = %0.2g, meanDwellTime = %0.2g',find(obj.sortMap == obj.Nspot),obj.meanTime2FirstEvent,obj.meanDwellTime),'Parent', obj.hras)
        end
        
          function plotTimePoint(obj,src,~)
            xcoords = get(obj.timeLine,'xdata');
            ii = round(xcoords(1));

            if size(obj.seq,3) > 3
                b = sqrt(obj.seq(obj.Nspot,:,4).^2+obj.seq(obj.Nspot,:,5).^2);
                title(sprintf('spot ID = %d, Frame = %d, \\Delta X = %0.2g, I = %0.2g, bg=%0.2g',obj.Nspot,ii,b(ii),obj.photons1I(obj.Nspot,ii),obj.bg1I(obj.Nspot,ii,2)*(obj.spotDetectorObject.paramsFit.PSFSigma*2)^2),'Parent',obj.hh);
            end              
            
            frames2 = obj.rawFitResultsCam2.Frame(obj.maskFilt2,:);
            frames1 = obj.rawFitResultsCam1.Frame(obj.maskFilt1,:);
    
            idx = find(frames2==ii);
            A(:,:,1+obj.flip) = squeeze(obj.rawFitResultsCam2.ROIStack(:,:,idx(obj.Nspot)));
      
            idx1 = find(frames1==ii);
            A(:,:,2-obj.flip) = squeeze(obj.rawFitResultsCam1.ROIStack(:,:,idx1(obj.Nspot)));
            xxyy = obj.rawFitResultsCam2.Coord(idx(obj.Nspot),:)-obj.rawFitResultsCam2.RoiStart(idx(obj.Nspot),:);
            xxyytarget = obj.coordsCam2(idx(obj.Nspot),1:2)-obj.rawFitResultsCam2.RoiStart(idx(obj.Nspot),:);

            gimr = mat2gray(A(:,:,1)',[min(obj.amin1,obj.amin2) max(obj.amax1,obj.amax2)]); %flip(flip(A(:,:,ii+1),1),2));
            
            if ~isempty(obj.himg1) && isvalid(obj.himg1)
                obj.himg1.CData=gimr*256;
            else
                obj.h1 = subplot(2,2,3,'Parent', obj.htab1);
                obj.himg1 = subimage(gimr);
                axis off;
                hold on
                title('Component','parent',obj.h1)
            end

            gimg = mat2gray(flip(flip(A(:,:,2),1),2),[min(obj.amin1,obj.amin2) max(obj.amax1,obj.amax2)]); 
            gimg = uint8(gimg * 256);
            gimr = uint8(gimr * 256);
            filler = zeros(size(gimr),'uint8');
            rgbImage = cat(3,gimr,gimg,filler);

            

            if ~isempty(obj.himg2) && isvalid(obj.himg2)
                obj.himg2.CData=rgbImage;
            else
                obj.h2 = subplot(2,2,4,'Parent', obj.htab1);
                obj.himg2 = subimage(rgbImage);
                axis off
                xlim([1 size(rgbImage,1)])
                ylim([1 size(rgbImage,2)])
                legend('Position')
                title('Reference & Component','parent',obj.h2)
            end
            drawnow
            hold on
            if length(obj.hp) >= 1 && ishandle(obj.hp(1))
                delete(obj.hp(1))
            end
            obj.hp(1) = plot(-0.1,-0.1,'h','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,'parent',obj.h2)
           
            if obj.colocalized
                if length(obj.hp) >= 2 && ishandle(obj.hp(2))
                    delete(obj.hp(2))
                end
            obj.hp(2) = plot(xxyy(1)+1,xxyy(2)+1,'h','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',8,'parent',obj.h2)
            end        
            if length(obj.hp) >= 3 && ishandle(obj.hp(3))
                delete(obj.hp(3))
            end
            obj.hp(3) = plot(xxyytarget(1)+1,xxyytarget(2)+1,'h','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',8,'parent',obj.h2)

            if obj.colocalized
                if length(obj.hp) >= 4 &&ishandle(obj.hp(4))
                    delete(obj.hp(4))
                end
            obj.hp(4) = plot(xxyy(1)+1,xxyy(2)+1,'h','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',8,'parent',obj.h1)
            end      
            if length(obj.hp) >= 5 &&ishandle(obj.hp(5))
                delete(obj.hp(5))
            end
            obj.hp(5) = plot(xxyytarget(1)+1,xxyytarget(2)+1,'h','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',8,'parent',obj.h1)
          end
                  
        function clickhline(obj,src,ev)
            clicked=get(gca,'currentpoint');
            obj.initialClick=clicked(1,2,1);
           
            set(gcf,'windowbuttonmotionfcn',@obj.draghline);
            set(gcf,'windowbuttonupfcn',@obj.dragdonehline);
        end
        
        function draghline(obj,src,ev)
          
            clicked=get(gca,'currentpoint');
            ycoord=clicked(1,2,1);
            set(obj.spotLine,'ydata',[ycoord ycoord]);
            if size(obj.seq,3) > 3
                a = obj.spotNrs(obj.sortMap);
                 obj.Nspot = a(round(ycoord));
                title(sprintf('spot ID = %d, meanTime2FirstEvent = %0.2g, meanDwellTime = %0.2g',obj.Nspot,obj.meanTime2FirstEvent,obj.meanDwellTime),'Parent', obj.hras)
            end            
        end

        function dragdonehline(obj,src,ev)
            clicked=get(gca,'currentpoint');
            finalClick=clicked(1,2,1);
            
            if finalClick == obj.initialClick
                ycoord = get(obj.spotLine,'ydata');
                a = obj.spotNrs(obj.sortMap);
                

                if finalClick < ycoord(1)
                    obj.Nspot = a(round(round(ycoord(1))-1));
                    set(obj.spotLine,'ydata',ycoord-1);
                else
                    obj.Nspot = a(round(round(ycoord(1))+1));
                    set(obj.spotLine,'ydata',ycoord+1);
                end
                title(sprintf('spot ID = %d, meanTime2FirstEvent = %0.2g, meanDwellTime = %0.2g',obj.Nspot,obj.meanTime2FirstEvent,obj.meanDwellTime),'Parent', obj.hras)
            end

            set(gcf,'windowbuttonmotionfcn','');
            set(gcf,'windowbuttonupfcn','');
            notify(obj,'spotChanged');
        end
    end
end