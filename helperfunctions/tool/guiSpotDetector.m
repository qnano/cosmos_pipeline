classdef guiSpotDetector < handle
   properties
    Figure
    handels
    Axes
    CLimit
    dataObject
    hTopTab
    htab1
    htab2
    lh 
    stringList
    coords
    cameraIndex
    NFramesDetect
    beginFrame
    endFrame
    addfilters
    him3d2
    him3d1
    hAxis3
    spotsClass
    
    htabgroup
    canCoordsCam1
     
    paramsPreFilterFits
    paramsFit
    paramsFilterFits      
    maskFilt1
    maskPreFiltCam1
    rawInitialFitResultsCam1
    detParCam1
    img1
    img2
    htab3
    
    hdragline
    
    hsub1 
    hsub2 
    hsub3 
    hsub4 
    hsub5
    hsub6
    l
    hAxis2
    
    darkROIs
    darkCoords
    
    parentMenu
    status=0;
    class = 'guiSpotDetector'
   end
   
    events
       ParamChange
       filterChange
   end
   
   methods
       function deletetry(obj,handle)
            try
                delete(handle)
            catch
            end
       end
       
       function whoAmI(obj,src,~)
            basevars = evalin('base','whos');
            testClassvars = basevars(strcmp({basevars.class},class(obj)));

            for i = 1:length(testClassvars)
                if(eq(evalin('base',testClassvars(i).name),obj))
                    obj.name =testClassvars(i).name;
                end
            end
       end
        function  settingspff(obj,src,~)
            prompt={'maximum circularity','minimum  circularity','maximum PH1','minimum PH1','pixel distance','maximum clustersize','minimum clustersize'};
            name = 'params Pre-FilterFits';
            defaultans = {num2str(obj.paramsPreFilterFits.circularityMax), num2str(obj.paramsPreFilterFits.circularityMin),...
                num2str(obj.paramsPreFilterFits.PH1Max), num2str(obj.paramsPreFilterFits.PH1Min),...
                 num2str(obj.paramsPreFilterFits.minPixelDist),num2str(obj.paramsPreFilterFits.clusterSizeMax), num2str(obj.paramsPreFilterFits.clusterSizeMin),...
                  ...
                  };
            options.Interpreter = 'tex';
            
            [answer] =  inputdlg(prompt,name,[1 50],defaultans,options);
            cancel = isempty(answer);
            if ~cancel
                obj.paramsPreFilterFits.circularityMax=str2num(answer{1});
                obj.paramsPreFilterFits.circularityMin=str2num(answer{2});
                obj.paramsPreFilterFits.PH1Max=str2num(answer{3});
                obj.paramsPreFilterFits.PH1Min= str2num(answer{4});
                obj.paramsPreFilterFits.minPixelDist = str2num(answer{5});
                obj.paramsPreFilterFits.clusterSizeMax = str2num(answer{6});
                obj.paramsPreFilterFits.clusterSizeMin = str2num(answer{7});
            end
        end
         function  settingsff(obj,src,~)
            prompt={'Iterations','MaxCudaFits','PSFSigma','BoxSize'};
            name = 'Fit params';
            defaultans = {num2str(obj.paramsFit.Iterations), num2str(obj.paramsFit.MaxCudaFits),...
                num2str(obj.paramsFit.PSFSigma), num2str(obj.paramsFit.BoxSize) };
            options.Interpreter = 'tex';
            
            [answer] =  inputdlg(prompt,name,[1 50],defaultans,options);
            cancel = isempty(answer);
            if ~cancel
                obj.paramsFit.Iterations=str2num(answer{1});
                obj.paramsFit.MaxCudaFits=str2num(answer{2});
                obj.paramsFit.PSFSigma= str2num(answer{3});
                obj.paramsFit.BoxSize = str2num(answer{4});
            end
         end
         
        
         function  settingsf(obj,src,~)
            prompt={'MinCRLBSTD','MinPValue','MinPhotons','MinBg','MaxCRLBSTD','MaxPValue','MaxPhotons','MaxBg','MinPixelDist'};
            name = 'Fit params';
            defaultans = {num2str(obj.paramsFilterFits.MinCRLBSTD), num2str(obj.paramsFilterFits.MinPValue),...
                 num2str(obj.paramsFilterFits.MinPhotons),...
            num2str(obj.paramsFilterFits.MinBg),...
            num2str(obj.paramsFilterFits.MaxCRLBSTD),...
            num2str(obj.paramsFilterFits.MaxPValue),...
            num2str(obj.paramsFilterFits.MaxPhotons),...
            num2str(obj.paramsFilterFits.MaxBg),...
            num2str(obj.paramsFilterFits.MinPixelDist) };
        
            options.Interpreter = 'tex';
            
            [answer] = inputdlg(prompt,name,[1 50],defaultans,options); 
            cancel = isempty(answer);
            if ~cancel
                obj.paramsFilterFits.MinCRLBSTD=str2num(answer{1});
                obj.paramsFilterFits.MinPValue=str2num(answer{2});
                obj.paramsFilterFits.MinPhotons=str2num(answer{3});
                obj.paramsFilterFits.MinBg=str2num(answer{4});
                obj.paramsFilterFits.MaxCRLBSTD=str2num(answer{5});
                obj.paramsFilterFits.MaxPValue=str2num(answer{6});
                obj.paramsFilterFits.MaxPhotons=str2num(answer{7});
                obj.paramsFilterFits.MaxBg=str2num(answer{8});
                obj.paramsFilterFits.MinPixelDist=str2num(answer{9});
            end
         end
        
      function obj = guiSpotDetector(h,dataObject,cameraIndex,NFramesDetect,name,nemuhandle,darkROIs,addfilters)
            if nargin < 3  || isempty(cameraIndex)
                cameraIndex = 1;
            end
            if nargin < 4  || isempty(NFramesDetect)
                NFramesDetect = 1;
            end
            if nargin < 7 || isempty(darkROIs)
                darkROIs=0;
            end
            if nargin < 8 || isempty(addfilters)
                obj.addfilters=true;
            else
                obj.addfilters=addfilters;
            end
            
            obj.status=0;
            
            obj.darkROIs=darkROIs;
            obj.NFramesDetect=NFramesDetect;
            obj.cameraIndex = cameraIndex;
            obj.dataObject = dataObject;

            if nargin < 5 || isempty(name)
                name = 'Spot Detector';                
            end
            obj.hTopTab = uitab(h, 'Title', name);
            
            if nargin > 5 && ~isempty(nemuhandle)
                obj.handels.parentMenu = uimenu('Parent',nemuhandle,'Label',name);
            else
                obj.handels.parentMenu = uimenu('Label',name);
            end

            uimenu(obj.handels.parentMenu,'Label','Fit Settings','Callback',@obj.settingsff);
            uimenu(obj.handels.parentMenu,'Label','Detect Spots','Separator','on','Callback',@obj.detectSpots);
            uimenu(obj.handels.parentMenu,'Label','Load','Separator','on','Callback',@obj.loadvars);
            uimenu(obj.handels.parentMenu,'Label','Save','Callback',@obj.savevars);
            if obj.addfilters
                obj.htabgroup = uitabgroup(obj.hTopTab);
                obj.htab1 = uitab(obj.htabgroup, 'Title', 'Detect Spots');
                obj.htab3 = uitab(obj.htabgroup, 'Title', 'Pre-Filter Spots');
                obj.htab2 = uitab(obj.htabgroup, 'Title', 'Filter Spots');
            else
                obj.htab1 = obj.hTopTab;
            end
            obj.paramsPreFilterFits = getDefaultParamsPreFilterFits;
            obj.paramsPreFilterFits.minPixelDist=7;
            obj.paramsPreFilterFits.clusterSizeMax=100;
            obj.paramsPreFilterFits.circularityMax=2;

            obj.paramsFit = getDefaultParamsFit;
            obj.paramsFit.FitSigma=true;

            obj.paramsFilterFits = getDefaultParamsFilterFits;
            obj.paramsFilterFits.MinPixelDist=2;
            
            obj.Figure = getParentFigure(h);  
            if obj.addfilters
                obj.stringList{1} = 'paramsFilterFits.MinPhotons';
                obj.stringList{2} = 'paramsFilterFits.MaxPhotons';
                obj.stringList{3} = 'paramsFilterFits.MinBg';
                obj.stringList{4} = 'paramsFilterFits.MaxBg';
                obj.stringList{5} = 'paramsFilterFits.MinPValue';
                obj.stringList{6} = 'paramsFilterFits.MaxPValue';
                obj.stringList{7} = 'paramsFilterFits.MinPixelDist';

                obj.stringList{8} = 'paramsPreFilterFits.circularityMin';
                obj.stringList{9} = 'paramsPreFilterFits.circularityMax';
                obj.stringList{10} = 'paramsPreFilterFits.clusterSizeMin';
                obj.stringList{11} = 'paramsPreFilterFits.clusterSizeMax';
                obj.stringList{12} = 'paramsPreFilterFits.PH1Min';
                obj.stringList{13} = 'paramsPreFilterFits.PH1Max';
                obj.stringList{14} = 'paramsPreFilterFits.minPixelDist';
            else
                obj.lh = [];
            end
            
            
      end
  
       function dataChange(obj,src,~)
           notify(obj,'ParamChange');
       end
        
        function vs = savevars(obj,src,~)
            vs.class = obj.class;
            vs.rawInitialFitResultsCam1 = obj.rawInitialFitResultsCam1;
            vs.maskPreFiltCam1 = obj.maskPreFiltCam1; 
            vs.maskFilt1 =  obj.maskFilt1; 
            
            vs.darkCoords = obj.darkCoords;
            vs.darkROIs = obj.darkROIs;
            vs.NFramesDetect = obj.NFramesDetect;
            
            vs.paramsPreFilterFits = obj.paramsPreFilterFits;
            vs.paramsFit = obj.paramsFit;
            vs.paramsFilterFits = obj.paramsFilterFits;
            vs.detParCam1 = obj.detParCam1;
            vs.canCoordsCam1 = obj.canCoordsCam1;
            
            if nargin > 1 && isobject(src)
                [file, pathname] = uiputfile([obj.class '.mat'] );
                if file ~= 0
                    answer{1} = fullfile(pathname, file);
                else
                    answer=[];
                end
                if ~isempty(answer)
                    save(answer{1},'vs','-v7.3')
                end
            end

        end

        function loadvars(obj,src,updateplot)
            
             if ishandle(src)
                [file, pathname] = uigetfile([obj.class '.mat'] );
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
                if nargin < 3 || isempty(updateplot) || isobject(updateplot)
                    updateplot=true;
                end
                updateplot=updateplot&(isempty(obj.dataObject) || obj.dataObject.status);

                if ischar(answer{1})
                    ws = load(answer{1});
                    ws=ws.vs;
                elseif isstruct(src)
                    ws=src;
                end

                class = obj.class;
                if ~strcmpi(ws.class,class)
                    error(['This is no ' obj.class ' class'])
                else           
                    obj.rawInitialFitResultsCam1 = ws.rawInitialFitResultsCam1;
                    obj.maskPreFiltCam1 = ws.maskPreFiltCam1;
                    obj.maskFilt1 =  ws.maskFilt1; 
                    obj.darkCoords = ws.darkCoords;
                    obj.darkROIs = ws.darkROIs;
                    obj.detParCam1 = ws.detParCam1;
                    obj.canCoordsCam1 = ws.canCoordsCam1;
                    obj.NFramesDetect = ws.NFramesDetect;

                    obj.paramsPreFilterFits = ws.paramsPreFilterFits;
                    obj.paramsFit = ws.paramsFit;
                    obj.paramsFilterFits = ws.paramsFilterFits;

                    if updateplot
                        obj.plotSpots;
                        obj.plotDarkSpots;                
                        obj.plotHistogramsFilter;
                        obj.plotHistogramsPreFilter;
                        obj.setFilterLines;
                    end
                    obj.status=1;
                end
            end
        end
       

        function  detectSpots(obj,src,~)
            if ~obj.dataObject.status
                error('Data for detection not loaded!')
            else

                if (isempty(obj.beginFrame) || isempty(obj.endFrame)) & ~isempty(obj.NFramesDetect)
                    obj.beginFrame = 1;
                    obj.endFrame = obj.NFramesDetect;
                elseif ~(~isempty(obj.beginFrame) && ~isempty(obj.endFrame))
                    error('which frame do you want to use for detection?')
                end
                [obj.canCoordsCam1,obj.detParCam1,~] = LLRMapv2(obj.dataObject.getData(obj.cameraIndex,[],[],obj.beginFrame:obj.endFrame) ,obj.paramsFit.PSFSigma); %dataCam{}(:,:,)
                obj.canCoordsCam1(:,3) = obj.canCoordsCam1(:,3)+obj.beginFrame-1;
                
                if obj.addfilters
                    [ obj.maskPreFiltCam1 ] =  preFilterFits(obj.canCoordsCam1,obj.detParCam1,obj.paramsPreFilterFits);
                else
                     obj.maskPreFiltCam1 = true(size(obj.canCoordsCam1,1),1);
                end
                

                obj.coords=round([obj.canCoordsCam1(:,2) obj.canCoordsCam1(:,1) obj.canCoordsCam1(:,3)]+(1.5*(2*obj.paramsFit.PSFSigma+1)-0.5).*[ones(size(obj.canCoordsCam1,1),1) ones(size(obj.canCoordsCam1,1),1) zeros(size(obj.canCoordsCam1,1),1)]);
    
                obj.fitSpots;
                obj.status=1;
            end
            
        end
      
       function plotSpots(obj,src,~)
         if ~isempty(obj.rawInitialFitResultsCam1)
            set(0,'CurrentFigure',obj.Figure);
            x1=0;
            y1=0;
            if isempty(obj.hAxis2)
                obj.hAxis2 =subplot(1,2,1,'Parent', obj.htab1);
                
                obj.him3d1= imtool3D(obj.hAxis2,double(stretch(obj.dataObject.getData(obj.cameraIndex,[],[],obj.beginFrame:obj.endFrame)))); %.dataCam{1}(:,:,))));
                obj.him3d1.pStretch.low = 0.01;
                obj.him3d1.pStretch.high = 0.999;
                hold on
            end
            for i=1:length(obj.spotsClass)
                obj.deletetry(obj.spotsClass(i))
            end
            
            obj.spotsClass(1) = plot(obj.rawInitialFitResultsCam1.Coord(~(obj.maskFilt1&obj.maskPreFiltCam1),1)+1,obj.rawInitialFitResultsCam1.Coord(~(obj.maskFilt1&obj.maskPreFiltCam1),2)+1,'xr','parent',obj.hAxis2);
            obj.spotsClass(2) = plot(obj.rawInitialFitResultsCam1.Coord(obj.maskFilt1&obj.maskPreFiltCam1,1)+1,obj.rawInitialFitResultsCam1.Coord(obj.maskFilt1&obj.maskPreFiltCam1,2)+1,'xg','parent',obj.hAxis2);

            if isempty(obj.hAxis3)
                obj.hAxis3 =subplot(1,2,2,'Parent', obj.htab1);
                
                obj.him3d2= imtool3D(obj.hAxis3,double(stretch(obj.dataObject.getData(obj.cameraIndex,[],[],obj.beginFrame:obj.endFrame)))); %.dataCam{1}(:,:,))));
                obj.him3d2.pStretch.low = 0.01;
                obj.him3d2.pStretch.high = 0.999;
            end
            
            if obj.darkROIs > 0
                DetectedCoords = [obj.rawInitialFitResultsCam1.Coord(:,1) obj.rawInitialFitResultsCam1.Coord(:,2)];

                gridcell=10; % This is the size of the dark AOI grid cell (smaller values lead to 

                % higher grid density of dark AOIs; this should be greater or equal to the AOI size
                mindist=10; % This is the minimal distance in pixels from non-dark AOIs to maintain 

                xmin=min(DetectedCoords(:,1))+gridcell/2;
                xmax=max(DetectedCoords(:,1))-gridcell/2;
                ymin=min(DetectedCoords(:,2))+gridcell/2;
                ymax=max(DetectedCoords(:,2))-gridcell/2;
                xgrid=linspace(xmin,xmax,(xmax-xmin)/gridcell);
                ygrid=linspace(ymin,ymax,(ymax-ymin)/gridcell);

                [x,y]=meshgrid(xgrid,ygrid);
                a=ones(size(x));

                for i=1:length(DetectedCoords(:,1)) 
                    xnotdark=find(abs(DetectedCoords(i,1)-xgrid)<mindist);
                    ynotdark=find(abs(DetectedCoords(i,2)-ygrid)<mindist);
                    a(ynotdark,xnotdark)=0;
                end
                obj.darkCoords=[];
                
                % Draw dark locations (blue) and RNA target locations (red)
                obj.darkCoords(:,1)=x(a==1);
                obj.darkCoords(:,2)=y(a==1);
                if obj.darkROIs > 1
                    idx = randperm(size(obj.darkCoords,1),min(size(obj.darkCoords,1),min(size(DetectedCoords,1),obj.darkROIs)));
                else
                    idx = randperm(size(obj.darkCoords,1),min(size(obj.darkCoords,1),size(DetectedCoords,1)));
                end
                obj.darkCoords = obj.darkCoords(idx,:);
                obj.spotsClass(3) = obj.plotDarkSpots;    
            end
         end         
       end
       
       function h = plotDarkSpots(obj,src,~)
           if ~isempty(obj.darkCoords)
               h = scatter(obj.darkCoords(:,1),obj.darkCoords(:,2), 'oy','Parent', obj.hAxis2)
           else
               h=[];
           end
       end
        
            function plotHistogramsPreFilter(obj,src,~)
                
              if ~isempty(obj.rawInitialFitResultsCam1)
 
                Nbin = floor(max(2,size(obj.rawInitialFitResultsCam1.Coord,1)/2));
                if Nbin>1
                    set(0,'CurrentFigure',obj.Figure);

                    obj.hsub2 =subplot(2,2,1,'Parent', obj.htab3);
                    histogram(obj.detParCam1.circularity,Nbin,'Normalization','probability','buttondownfcn',@obj.clickline);
                    title('Circularity');

                    obj.hsub3 =subplot(2,2,2,'Parent', obj.htab3);                
                    histogram(obj.detParCam1.clusterSize,Nbin,'Normalization','probability','buttondownfcn',@obj.clickline);
                    title('Cluster Size');
                    

                    obj.hsub4 =subplot(2,2,3,'Parent', obj.htab3);                
                    histogram(obj.detParCam1.PH1,Nbin,'Normalization','probability','buttondownfcn',@obj.clickline);
                    title('P_{H1}');                  
                    
                    obj.hsub6 =subplot(2,2,4,'Parent', obj.htab3,'buttondownfcn',@obj.clickline);
                    [~,dist] = knnsearch(obj.canCoordsCam1, obj.canCoordsCam1,'k',2);   %% for multiple frames multiply coords by .*repmat([1 1 2*obj.paramsFilterFits.MinPixelDist],[size(obj.rawInitialFitResultsCam1.Coord,1) 1])
                    histogram(dist(:,2),Nbin,'Normalization','probability','buttondownfcn',@obj.clickline);
                    title('\Delta x Neighbour');
                    
                    drawnow;
                    obj.l(8)=line([1 1].*obj.paramsPreFilterFits.circularityMin,get(obj.hsub2, 'YLim'),'Color','g','Linewidth',2,'parent',obj.hsub2,'tag','line1','hittest','off'); %min(obj.rawInitialFitResultsCam1.Photons)
                    drawnow;
                    obj.l(9)=line([1 1].*obj.paramsPreFilterFits.circularityMax,get(obj.hsub2, 'YLim'),'Color','r','Linewidth',2,'parent',obj.hsub2,'tag','line2','hittest','off');
                    
                    obj.l(10)=line([1 1].*obj.paramsPreFilterFits.clusterSizeMin,get(obj.hsub3, 'YLim'),'Color','g','Linewidth',2,'parent',obj.hsub3,'tag','line1','hittest','off'); %min(obj.rawInitialFitResultsCam1.Bg)
                    obj.l(11)=line([1 1].*obj.paramsPreFilterFits.clusterSizeMax,get(obj.hsub3, 'YLim'),'Color','r','Linewidth',2,'parent',obj.hsub3,'tag','line2','hittest','off');
                    
                    obj.l(12)=line([1 1].*obj.paramsPreFilterFits.PH1Min,get(obj.hsub4, 'YLim'),'Color','g','Linewidth',2,'parent',obj.hsub4,'tag','line1','hittest','off');
                    obj.l(13)=line([1 1].*obj.paramsPreFilterFits.PH1Max,get(obj.hsub4, 'YLim'),'Color','r','Linewidth',2,'parent',obj.hsub4,'tag','line2','hittest','off'); %max(PFA_adj)

                    obj.l(14)=line([1 1].*obj.paramsPreFilterFits.minPixelDist,get(obj.hsub6, 'YLim'),'Color','g','Linewidth',2,'parent',obj.hsub6,'tag','line1','hittest','off');%min(dist(:,2))
                    
                    set(obj.hsub1 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');                    
                    set(obj.hsub2 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');                    
                    set(obj.hsub3 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');                    
                    set(obj.hsub4 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');
                    set(obj.hsub5 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');                   
                    set(obj.hsub6 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');
                else
                    delete(obj.hsub1);
                    delete(obj.hsub2);
                    delete(obj.hsub3);
                    delete(obj.hsub4);
                end
             end    
         end
          
         function plotHistogramsFilter(obj,src,~)
              if ~isempty(obj.rawInitialFitResultsCam1)
 
                Nbin = floor(max(2,size(obj.rawInitialFitResultsCam1.Coord,1)/2));
                if Nbin>1
                    set(0,'CurrentFigure',obj.Figure);
                    
                    obj.hsub2 =subplot(2,2,1,'Parent', obj.htab2);
                    histogram(obj.rawInitialFitResultsCam1.Photons,Nbin,'Normalization','probability','buttondownfcn',@obj.clickline);
                    title('Intensity');

                    obj.hsub3 =subplot(2,2,2,'Parent', obj.htab2);                
                    histogram(obj.rawInitialFitResultsCam1.Bg,Nbin,'Normalization','probability','buttondownfcn',@obj.clickline);
                    title('Background');
                    
                    [ ~,~,PFA_adj ]=fdr_bh(reshape(obj.rawInitialFitResultsCam1.LL(3,:),...
                    [prod(size(obj.rawInitialFitResultsCam1.LL(3,:))) 1]),min(0.05,obj.paramsFilterFits.MaxPValue),'dep','no');
                    PFA_adj=min(PFA_adj,1);
                    obj.hsub4 =subplot(2,2,3,'Parent', obj.htab2,'buttondownfcn',@obj.clickline);
                    histogram(PFA_adj,Nbin,'Normalization','probability','buttondownfcn',@obj.clickline);
                    title('P_{FA}');
                   
                    drawnow;
                    obj.hsub6 =subplot(2,2,4,'Parent', obj.htab2,'buttondownfcn',@obj.clickline);
                    [~,dist] = knnsearch(obj.rawInitialFitResultsCam1.Coord, obj.rawInitialFitResultsCam1.Coord,'k',2);   %% for multiple frames multiply coords by .*repmat([1 1 2*obj.paramsFilterFits.MinPixelDist],[size(obj.rawInitialFitResultsCam1.Coord,1) 1])
                    histogram(dist(:,2),Nbin,'Normalization','probability','buttondownfcn',@obj.clickline);
                    title('\Delta x Neighbour');
                    
                    obj.l(2)=line([1 1].*min(obj.paramsFilterFits.MaxPhotons,max(obj.rawInitialFitResultsCam1.Photons)),get(obj.hsub2, 'YLim'),'Color','r','Linewidth',2,'parent',obj.hsub2,'tag','line2','hittest','off');
                    obj.l(1)=line([1 1].*obj.paramsFilterFits.MinPhotons,get(obj.hsub2, 'YLim'),'Color','g','Linewidth',2,'parent',obj.hsub2,'tag','line1','hittest','off'); %min(obj.rawInitialFitResultsCam1.Photons)
                    delete(obj.l(2))
                    obj.l(2)=line([1 1].*min(obj.paramsFilterFits.MaxPhotons,max(obj.rawInitialFitResultsCam1.Photons)),get(obj.hsub2, 'YLim'),'Color','r','Linewidth',2,'parent',obj.hsub2,'tag','line2','hittest','off');
                    
                    obj.l(3)=line([1 1].*obj.paramsFilterFits.MinBg,get(obj.hsub3, 'YLim'),'Color','g','Linewidth',2,'parent',obj.hsub3,'tag','line1','hittest','off'); %min(obj.rawInitialFitResultsCam1.Bg)
                    obj.l(4)=line([1 1].*min(obj.paramsFilterFits.MaxBg,max(obj.rawInitialFitResultsCam1.Bg)),get(obj.hsub3, 'YLim'),'Color','r','Linewidth',2,'parent',obj.hsub3,'tag','line2','hittest','off');
                    
                    obj.l(5)=line([1 1].*obj.paramsFilterFits.MinPValue,get(obj.hsub4, 'YLim'),'Color','g','Linewidth',2,'parent',obj.hsub4,'tag','line1','hittest','off');
                    obj.l(6)=line([1 1].*obj.paramsFilterFits.MaxPValue,get(obj.hsub4, 'YLim'),'Color','r','Linewidth',2,'parent',obj.hsub4,'tag','line2','hittest','off'); %max(PFA_adj)

                    obj.l(7)=line([1 1].*obj.paramsFilterFits.MinPixelDist,get(obj.hsub6, 'YLim'),'Color','g','Linewidth',2,'parent',obj.hsub6,'tag','line1','hittest','off');%min(dist(:,2))
                    
                    set(obj.hsub1 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');                    
                    set(obj.hsub2 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');                    
                    set(obj.hsub3 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');                    
                    set(obj.hsub4 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');
                    set(obj.hsub5 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');                   
                    set(obj.hsub6 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');
                else
                    delete(obj.hsub1);
                    delete(obj.hsub2);
                    delete(obj.hsub3);
                    delete(obj.hsub4);
                end
             end    
         end
         
          function clickline(obj,src,ev)
            clicked=get(gca,'currentpoint');
            xcoord=clicked(1,1,1);
            obj.hdragline.l1=findobj(get(gca,'Children'),'tag','line1');
            if ~isempty(obj.hdragline.l1)
                xcoord1 = get(obj.hdragline.l1,'xdata');
            else
                xcoord1=Inf.*[1 1];
            end
            obj.hdragline.l2=findobj(get(gca,'Children'),'tag','line2');
            if ~isempty(obj.hdragline.l2)
                xcoord2 = get(obj.hdragline.l2,'xdata');
            else
                xcoord2=Inf.*[1 1];
            end
            if sum(abs(xcoord-xcoord1)) <= sum(abs(xcoord-xcoord2))
                obj.hdragline.current = obj.hdragline.l1;
            else
                obj.hdragline.current=obj.hdragline.l2;
            end
            
            set(gcf,'windowbuttonmotionfcn',@obj.dragline)
            set(gcf,'windowbuttonupfcn',@obj.dragdone)
          end
        
        function dragline(obj,src,ev)
            clicked=get(gca,'currentpoint');
            xcoord=clicked(1,1,1);
            set(obj.hdragline.current,'xdata',[xcoord xcoord]);
            if all(~isempty(obj.hdragline.l2)) && all(~isempty(obj.hdragline.l1))
                if (get(obj.hdragline.l2,'xdata') < get(obj.hdragline.l1,'xdata'))
                 set(obj.hdragline.l1,'xdata',[xcoord xcoord]);
                 set(obj.hdragline.l2,'xdata',[xcoord xcoord]);
                end
            end
        end


        function dragdone(obj,src,ev)
            set(gcf,'windowbuttonmotionfcn','');
            set(gcf,'windowbuttonupfcn','');
            obj.hdragline=[];
            obj.filterUpdate;
        end

        function  fitSpots(obj,src,~)
            if ~isempty(obj.dataObject.getData(obj.cameraIndex))
                obj.deletetry(obj.hAxis2);
                obj.deletetry(obj.hAxis3);
                obj.hAxis2 = [];
                obj.hAxis3 = [];
                
                obj.paramsFit.FitSigma=true;

                [ obj.rawInitialFitResultsCam1 ] = fitBoxCenters( single(squeeze(obj.dataObject.getData(obj.cameraIndex))),[obj.coords],obj.paramsFit);
                
                %some logic to kick out beads or other feducial markers
                obj.paramsFilterFits.MaxPhotons = max(obj.rawInitialFitResultsCam1.Photons);%...
                obj.filter;
                obj.plotSpots;
                if obj.addfilters
                    obj.plotHistogramsFilter;
                    obj.plotHistogramsPreFilter;
                end
            end
        end
        
        function filter(obj,src,~)  
            if obj.addfilters
                [ obj.maskPreFiltCam1 ] =  preFilterFits(obj.canCoordsCam1,obj.detParCam1,obj.paramsPreFilterFits);
                [ obj.maskFilt1 ] =  filterFits(obj.rawInitialFitResultsCam1,obj.paramsFilterFits);
            else
                 obj.maskPreFiltCam1 = true(size(obj.canCoordsCam1,1),1);
                 obj.maskFilt1 = true(size(obj.rawInitialFitResultsCam1,1),1);
            end
        end
        
        function filterUpdate(obj,src,~)
            for i=1:length(obj.l)
                xdata = get(obj.l(i),'xdata');
                eval(['obj.' obj.stringList{i} ' = ' num2str(xdata(1)) ';']);
             end
            obj.filter;
            obj.plotSpots;
            notify(obj,'filterChange');
        end
        
        function setFilterLines(obj,src,~)
            for i=1:length(obj.l)
                eval(['set(obj.l(i),''xdata'',[1 1]*obj.' obj.stringList{i} ');']);
             end
        end
        
        function rs = getSummary(obj)
              
            if isempty(obj.disttFormEst)
                rs = {'Estimate alignment first'};
            else
                 txt1 = sprintf('Detectedt %d spots\n',sum(obj.maskFilt1 & obj.maskPreFiltCam1));
                rs = {txt1};
            end
        end
        
        function clear(obj,src,~)
            obj.paramsPreFilterFits=[];
            obj.paramsFit=[];
            obj.paramsFilterFits=[];
              
            delete(obj.img1)
            delete(obj.img2)
        end
   end
end