classdef spotDetectorMultiplexRNAWrapper < handle
   properties
    Figure
    him3d
    spotDetectorObjects
    maskPreFiltCam1
    maskFilt1
    lh2
    lh
    hsub1
    hsub2
    hsub5
    hsub6
    hsub4
    hsub3
    stringList
    htab1
    htab3
    htab2
    points
    hdragline
    NFramesDetect
    pointsExcluded
    addfilters
    C
    l
    htabgroup
    driftObject
    detParCam1=[];
    canCoordsCam1=[];
    rawInitialFitResultsCam1=[];    
    paramsPreFilterFits
    paramsFilterFits
    paramsFit
    xdata

    darkCoords
    darkROIs
    newCoords
    dataObject
    status
    class = 'spotDetectorMultiplexRNAWrapper';
   end
   
    events
       filterChange
       ParamChange
       exclusionChange
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
       function setFrameNumbers(obj,C,~)
           answer{1}='1';
           answer{2}='2';
           if size(C,2) == 2 && size(C,1) == length(obj.spotDetectorObjects)
                obj.C=C;
           else
                for i=1:length(obj.spotDetectorObjects)
                    prompt={['begin frame RNA ' num2str(i)],['end frame RNA ' num2str(i)]};
                    name = ['detection frames RNA ' num2str(i)];
                    defaultans = { answer{1}, answer{2} };
                    options.Interpreter = 'tex';

                    [answer] =  inputdlg(prompt,name,[1 50],defaultans,options);
                    cancel = isempty(answer);
                    if ~cancel
                        C(i,1) =str2num(answer{1});
                        C(i,2) =str2num(answer{2});
                        obj.C=C;
                    else
                        error('This was not supposed to happen');
                    end                
                end
           end
           if ~isempty(obj.C)
               for i=1:length(obj.spotDetectorObjects)
                    obj.spotDetectorObjects{i}.beginFrame = C(i,1);
                    obj.spotDetectorObjects{i}.endFrame = C(i,2);
               end
           end 
       end

         
  
      function obj = spotDetectorMultiplexRNAWrapper(h,gdlObj,driftObject,MultiplexRNA,cameraIndex)
          obj.addfilters = false; %this is still buggy
          obj.darkROIs=true;
        if nargin < 5 || isempty(cameraIndex)
            cameraIndex=1;
        end
        addSpotDetectorFilters = true;
        obj.status=0;
        obj.Figure = getParentFigure(h); 
        obj.dataObject=gdlObj; 
        obj.driftObject=driftObject;
        topTab = uitab(h, 'Title', 'Spot Detector');
        htabgroup = uitabgroup(topTab);
        f = uimenu('Label','Spot Detector');
        for i=1:MultiplexRNA
            obj.spotDetectorObjects{i} = guiSpotDetector(htabgroup,gdlObj,cameraIndex,1,['RNA ' num2str(i)],f,[],addSpotDetectorFilters);
            obj.deletetry(obj.spotDetectorObjects{i}.handels.parentMenu)
        end

        if  obj.addfilters
            subTopTab = uitab(htabgroup, 'Title', 'Different Species');
            obj.htabgroup = uitabgroup(subTopTab);
            obj.htab1 = uitab(obj.htabgroup, 'Title', 'Detected Spots');
            obj.htab3 = uitab(obj.htabgroup, 'Title', 'Pre-Filter Spots');
            obj.htab2 = uitab(obj.htabgroup, 'Title', 'Filter Spots');
        else
            obj.htab1 = uitab(htabgroup, 'Title', 'Different Species');
        end
        uimenu(f,'Label','Fit Settings','Callback',@obj.settingsff);  
        uimenu(f,'Label','Set Detection Frames','Callback',@obj.setFrameNumbers);       
        uimenu(f,'Label','Detect Spots','Separator','on','Callback',@obj.detectSpots);
          
        uimenu(f,'Label','Load','Separator','on','Callback',@obj.loadvars);  
        uimenu(f,'Label','Save','Callback',@obj.savevars); 
        obj.paramsFit = getDefaultParamsFit;
        obj.paramsFit.FitSigma=true;
        
        if obj.addfilters
            obj.paramsPreFilterFits = getDefaultParamsPreFilterFits;
            obj.paramsPreFilterFits.minPixelDist=7;
            obj.paramsPreFilterFits.clusterSizeMax=100;
            obj.paramsPreFilterFits.circularityMax=2;

            obj.paramsFilterFits = getDefaultParamsFilterFits;
            obj.paramsFilterFits.MinPixelDist=2;

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
        end
        if addSpotDetectorFilters
            for i=1:length(obj.spotDetectorObjects)
                obj.lh2{i} = addlistener(obj.spotDetectorObjects{i},'filterChange',@obj.subFilterUpdate)
            end  
        end
        
      end
      function subFilterUpdate(obj,src,~)
        obj.status=0;
        obj.calculateExclusion([],true);
        for i=1:length(obj.spotDetectorObjects)
            obj.spotDetectorObjects{i}.darkCoords=obj.darkCoords;               
        end
        obj.status=1;
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
                for i=1:length(obj.spotDetectorObjects)
                    obj.spotDetectorObjects{i}.paramsFit.Iterations=str2num(answer{1});
                    obj.spotDetectorObjects{i}.paramsFit.MaxCudaFits=str2num(answer{2});
                    obj.spotDetectorObjects{i}.paramsFit.PSFSigma= str2num(answer{3});
                    obj.spotDetectorObjects{i}.paramsFit.BoxSize = str2num(answer{4});
                end
            end
      end
         
          
      
      function detectSpots(obj,src,~)
          if ~obj.driftObject.status || ~obj.dataObject.status
              error('No drift estimated or data loaded')
          end
              
        for i=1:length(obj.spotDetectorObjects)
            if ~isempty(obj.C)
                obj.spotDetectorObjects{i}.detectSpots(obj.C(i,:));
            else
                obj.setFrameNumbers([]);
                obj.spotDetectorObjects{i}.detectSpots(obj.C(i,:));
            end
        end
        obj.calculateExclusion;
        for i=1:length(obj.spotDetectorObjects)
            obj.spotDetectorObjects{i}.darkCoords=obj.darkCoords;
            obj.spotDetectorObjects{i}.plotDarkSpots;
        end
        if obj.addfilters
            obj.detParCam1.circularity=[];
            obj.detParCam1.clusterSize=[];
            obj.detParCam1.PH1=[];
            obj.canCoordsCam1=[];
            obj.rawInitialFitResultsCam1.Photons=[];
            obj.rawInitialFitResultsCam1.Bg=[];
            obj.rawInitialFitResultsCam1.LL=[];
            obj.rawInitialFitResultsCam1.Coord=[];
            obj.rawInitialFitResultsCam1.CRLB_STD= [];
            obj.rawInitialFitResultsCam1.Sigma = [];
            obj.rawInitialFitResultsCam1.RoiStart=[];
            obj.rawInitialFitResultsCam1.ROIStack=[];
            obj.rawInitialFitResultsCam1.Frame= [];
            for i=1:length(obj.spotDetectorObjects)
                obj.detParCam1.circularity = cat(1,obj.detParCam1.circularity,obj.spotDetectorObjects{i}.detParCam1.circularity); 
                obj.detParCam1.clusterSize = cat(1,obj.detParCam1.clusterSize,obj.spotDetectorObjects{i}.detParCam1.clusterSize); 
                obj.detParCam1.PH1 = cat(2,obj.detParCam1.PH1,obj.spotDetectorObjects{i}.detParCam1.PH1); 
                obj.canCoordsCam1 = cat(1,obj.canCoordsCam1,obj.spotDetectorObjects{i}.canCoordsCam1);  
                obj.rawInitialFitResultsCam1.Photons= cat(1,obj.rawInitialFitResultsCam1.Photons,obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Photons); 
                obj.rawInitialFitResultsCam1.Bg= cat(1,obj.rawInitialFitResultsCam1.Bg,obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Bg); 
                obj.rawInitialFitResultsCam1.LL= cat(2,obj.rawInitialFitResultsCam1.LL,obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.LL);
                drift = repmat(obj.driftObject.drift(i,:),size(obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Coord,1),1);
                obj.rawInitialFitResultsCam1.Coord = cat(1,obj.rawInitialFitResultsCam1.Coord,obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Coord+drift); 
                obj.rawInitialFitResultsCam1.CRLB_STD= cat(1,obj.rawInitialFitResultsCam1.CRLB_STD,obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.CRLB_STD); 
                obj.rawInitialFitResultsCam1.Sigma= cat(1,obj.rawInitialFitResultsCam1.Sigma,obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Sigma); 
                obj.rawInitialFitResultsCam1.RoiStart= cat(1,obj.rawInitialFitResultsCam1.RoiStart,obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.RoiStart+drift); 
                obj.rawInitialFitResultsCam1.ROIStack=cat(3,obj.rawInitialFitResultsCam1.ROIStack,obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.ROIStack); 
                obj.rawInitialFitResultsCam1.Frame= cat(1,obj.rawInitialFitResultsCam1.Frame,obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Frame); 
            end
            obj.plotHistogramsFilter;
            obj.plotHistogramsPreFilter;
            obj.setFilterLines;
        end
        obj.status=1;
      end
      
            function plotHistogramsPreFilter(obj,src,~)
                
              if ~isempty(obj.detParCam1)
 
                Nbin = floor(max(2,size(obj.canCoordsCam1,1)/2));
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
                    title('PH1');                  

                    obj.hsub6 =subplot(2,2,4,'Parent', obj.htab3,'buttondownfcn',@obj.clickline);
                    [~,dist] = knnsearch(obj.canCoordsCam1, obj.canCoordsCam1,'k',2);   %% for multiple frames multiply coords by .*repmat([1 1 2*obj.paramsFilterFits.MinPixelDist],[size(obj.detParCam1,1) 1])
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
                    
                    set(0,'CurrentFigure',obj.Figure);
                    obj.l(1)=line([1 1].*obj.paramsFilterFits.MinPhotons,get(obj.hsub2, 'YLim'),'Color','g','Linewidth',2,'parent',obj.hsub2,'tag','line1','hittest','off'); %min(obj.rawInitialFitResultsCam1.Photons)                    

                    obj.paramsFilterFits.MaxPhotons = max(obj.rawInitialFitResultsCam1.Photons);
                    obj.paramsFilterFits.MaxBg = max(obj.rawInitialFitResultsCam1.Bg);
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

      
      function setFilterLines(obj,src,~)
            for i=1:length(obj.l)
                eval(['set(obj.l(i),''xdata'',[1 1]*obj.' obj.stringList{i} ');']);
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
      
      function filterUpdate(obj,src,~)
        if obj.status == 1
            obj.status=0;
            for i=1:length(obj.l)
                xdata = get(obj.l(i),'xdata');
                eval(['obj.' obj.stringList{i} ' = ' num2str(xdata(1)) ';']);
            end
            [ obj.maskPreFiltCam1 ] =  preFilterFits(obj.canCoordsCam1,obj.detParCam1,obj.paramsPreFilterFits);
            [ obj.maskFilt1 ] =  filterFits(obj.rawInitialFitResultsCam1,obj.paramsFilterFits);
            idx=1;
             for i=1:length(obj.spotDetectorObjects)
                 obj.spotDetectorObjects{i}.maskPreFiltCam1 = obj.maskPreFiltCam1(idx:idx+size(obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Coord,1)-1);
                 obj.spotDetectorObjects{i}.maskFilt1 = obj.maskFilt1(idx:idx+size(obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Coord,1)-1);
                 idx=idx+size(obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Coord,1)-1;
             end
            obj.calculateExclusion([],true);
            for i=1:length(obj.spotDetectorObjects)
                obj.spotDetectorObjects{i}.darkCoords=obj.darkCoords;
            end
            obj.status=1;
        end
      end
          
      
      function calculateExclusion(obj,src,updateplot)  
          if nargin < 3 || isempty(updateplot) || isobject(updateplot)
                updateplot=true;
            end
            
          detectedPoints=[];
          thres=5;
          for i=1:length(obj.spotDetectorObjects)
             mask = obj.spotDetectorObjects{i}.maskPreFiltCam1 & obj.spotDetectorObjects{i}.maskFilt1;
             if ~isempty(detectedPoints) 
                drift = repmat(obj.driftObject.drift(i,:),size(obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Coord,1),1);
                obj.newCoords{i} = obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Coord+drift;
                [idx,dist] = knnsearch(detectedPoints,obj.newCoords{i});
                
                indexNewPoints = (dist > thres) & mask;
                obj.spotDetectorObjects{i}.maskPreFiltCam1 = indexNewPoints;
                obj.spotDetectorObjects{i}.maskFilt1 = indexNewPoints;
                obj.points{i} = obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Coord(indexNewPoints,:);
                obj.pointsExcluded{i} = obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Coord((dist > thres) & ~mask,:);
                detectedPoints = cat(1,detectedPoints,obj.newCoords{i}(:,:));
             else
                obj.newCoords{i}=obj.spotDetectorObjects{i}.rawInitialFitResultsCam1.Coord;
                detectedPoints = cat(1,detectedPoints,obj.newCoords{i}(mask,:));
                obj.points{i} = detectedPoints;
                obj.pointsExcluded{i} = obj.newCoords{i}(~mask,:);
            end
          end
        if obj.darkROIs
            gridcell=10; % This is the size of the dark AOI grid cell (smaller values lead to 

            % higher grid density of dark AOIs; this should be greater or equal to the AOI size
            mindist=10; % This is the minimal distance in pixels from non-dark AOIs to maintain 
            DetectedCoords=[];
            for i=1:length(obj.spotDetectorObjects)                   
                DetectedCoords = cat(1,DetectedCoords, obj.newCoords{i});
            end
            maxSize = size(obj.spotDetectorObjects{end}.dataObject.getData(1),1);

            xmin=max(min(DetectedCoords(:,1)),1)+gridcell/2;
            xmax=min(maxSize,max(DetectedCoords(:,1)))-gridcell/2;
            ymin=max(min(DetectedCoords(:,2)),1)+gridcell/2;
            ymax=min(max(DetectedCoords(:,2)),maxSize)-gridcell/2;
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
        end
          if updateplot
            obj.plotImages;
          end
      end
      
      
        function plotImages(obj,src,~)
            if ~isempty(obj.points) && ~ isempty(obj.spotDetectorObjects{end}.dataObject.getData(1))
                maxSize = size(obj.spotDetectorObjects{end}.dataObject.getData(1),1);
                subAxis =subplot(1,2,1,'Parent', obj.htab1);
                
                c=['c','b','g','m'];
                obj.deletetry(obj.him3d);
                obj.him3d= imtool3D(subAxis,double(stretch(obj.spotDetectorObjects{end}.dataObject.getData(1,[],[],obj.C(1,1):obj.C(end,2))))); %.dataCam{1}(:,:,))));
                obj.him3d.pStretch.low = 0.01;
                obj.him3d.pStretch.high = 0.999;
                obj.him3d.showSlice(round(get(obj.him3d.Slider,'value')));
                hold on
               
                for i=1:length(obj.spotDetectorObjects)                   
                    plot(obj.him3d.haxis,min(maxSize,max(obj.points{i}(:,1)+1,0)),min(maxSize,max(obj.points{i}(:,2)+1,0)),[c(i) 'x'] ) %,'marker'
                    hold on
                    legendtexts{i}=['RNA ' num2str(i)];
                end
                if obj.darkROIs
                    legendtexts{end+1} = 'Dark Regions';
                    scatter(min(maxSize,max(obj.darkCoords(:,1),0)),min(maxSize,max(obj.darkCoords(:,2),0)), 'oy','Parent', obj.him3d.haxis)                
                end
                 legendtexts{end+1} = 'Excluded Spots';
                for i=1:length(obj.spotDetectorObjects)
                    plot(obj.him3d.haxis,min(maxSize,max(obj.pointsExcluded{i}(:,1)+1,0)),min(maxSize,max(obj.pointsExcluded{i}(:,2)+1,0)),[ 'xr'] ) %,'marker'
                end
                
                legend(legendtexts);          
            end    
        end
        
        function vs = savevars(obj,src,~)
          
            vs.class = obj.class;
            vs.darkCoords = obj.darkCoords;
            vs.darkROIs = obj.darkROIs;
            vs.newCoords = obj.newCoords;
            vs.points = obj.points;
            vs.NFramesDetect = obj.NFramesDetect;
            vs.C = obj.C;
            
            for i=1:length(obj.spotDetectorObjects)
                vs.rawInitialFitResultsCam1{i} = obj.spotDetectorObjects{i}.rawInitialFitResultsCam1;
                vs.maskPreFiltCam1{i} = obj.spotDetectorObjects{i}.maskPreFiltCam1;
                vs.maskFilt1{i} = obj.spotDetectorObjects{i}.maskFilt1;
                vs.beginFrame{i} = obj.spotDetectorObjects{i}.beginFrame;
                vs.endFrame{i} = obj.spotDetectorObjects{i}.endFrame;
                vs.detParCam1{i} = obj.spotDetectorObjects{i}.detParCam1;
                vs.canCoordsCam1{i} = obj.spotDetectorObjects{i}.canCoordsCam1;
            end
            
            if nargin > 1 && isobject(src)
                prompt={'File name'};
                name = 'Save';
                defaultans = {['spotDetectorMultiplexRNAWrapper' date '.mat']};
                options.Interpreter = 'tex';
                answer = inputdlg(prompt,name,[1 40],defaultans,options);
                if ~isempty(answer)
                    save(answer{1},'vs','-v7.3')
                end
            end
        end

        function loadvars(obj,src,updateplot)
            
            if ishandle(src)
                prompt={'File name'};
                name = 'Save';
                defaultans = {['spotDetectorMultiplexRNAWrapper' date '.mat']};
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

                if ischar(answer{1})
                    ws = load(answer{1});
                    ws=ws.vs;
                elseif isstruct(src)
                    ws=src;
                end
                if ~strcmpi(ws.class,obj.class)
                    error(['This is no ' obj.class ' class'])
                else
                    obj.darkCoords = ws.darkCoords;
                    obj.darkROIs = ws.darkROIs;
                    obj.newCoords = ws.newCoords;
                    obj.points = ws.points;
                    obj.NFramesDetect = ws.NFramesDetect;
                    obj.C = ws.C;
                    for i=1:length(obj.spotDetectorObjects)
                        obj.spotDetectorObjects{i}.rawInitialFitResultsCam1 = ws.rawInitialFitResultsCam1{i} ;
                        obj.spotDetectorObjects{i}.maskPreFiltCam1 = ws.maskPreFiltCam1{i};
                        obj.spotDetectorObjects{i}.maskFilt1 =  ws.maskFilt1{i};
                        obj.spotDetectorObjects{i}.beginFrame = ws.beginFrame{i} ;
                        obj.spotDetectorObjects{i}.endFrame = ws.endFrame{i};
                        obj.spotDetectorObjects{i}.detParCam1 = ws.detParCam1{i};
                        obj.spotDetectorObjects{i}.canCoordsCam1 = ws.canCoordsCam1{i};
                    end


                    obj.calculateExclusion([],updateplot);
                    for i=1:length(obj.spotDetectorObjects)
                        obj.spotDetectorObjects{i}.darkCoords=obj.darkCoords;
                        if updateplot
                            obj.spotDetectorObjects{i}.plotSpots;
                            obj.spotDetectorObjects{i}.plotDarkSpots;
                            obj.spotDetectorObjects{i}.plotHistogramsFilter;
                            obj.spotDetectorObjects{i}.plotHistogramsFilter;
                            obj.spotDetectorObjects{i}.plotHistogramsPreFilter;
                            obj.spotDetectorObjects{i}.setFilterLines;
                        end

                    end
                end
                obj.status=1;
                for i=1:length(obj.spotDetectorObjects)
                    obj.spotDetectorObjects{i}.status=1;
                end
            end
        end
   end
end