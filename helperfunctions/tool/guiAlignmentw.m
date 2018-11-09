classdef guiAlignmentw < handle
   properties
    %handles
    Figure
    Axes
    CLimit
    dataObject
    htab1
    htab2
    topTab
    htabgroup
    htabgroup2
    
    %data
    spotDetectorCam1
    spotDetectorCam2
    
    %settings
    paramsPreFilterFits
    paramsFit
    paramsFilterFits      
    
    %results
    mGellerCam1
    mGellerCam2
    cpPoints
    tformNoised
    cpPointsMoved
    cpPointsUnMoved
    tformTotal
    disttFormEst
    initialtransVec
    idx
    img1
    img2
    rgbImage
    mseFormEst
    
    handels
    
    status = 0;
    class ='guiAlignmentw';
   end
   
    events
       ParamChange
   end
   
   methods
       
       function whoAmI(obj,src,~)
            basevars = evalin('base','whos');
            testClassvars = basevars(strcmp({basevars.class},class(obj)));

            found = false;
            for i = 1:length(testClassvars)
                if(eq(evalin('base',testClassvars(i).name),obj))
                    found = true;
                    obj.name =testClassvars(i).name;
                end
            end
       end
       
       function deletetry(obj,handle)
            try
                delete(handle)
            catch
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
                obj.paramsFit.clusterSizeMin = str2num(answer{5});
                obj.paramsFit.clusterSizeMax = str2num(answer{6});  
            end
         end
         
        
        
        
      function obj = guiAlignmentw(h,dataObject)
            obj.dataObject = dataObject;
            nFrames=3
            
            obj.topTab = uitab(h, 'Title', 'Alignment');
            obj.htabgroup = uitabgroup(obj.topTab);
            obj.htab1 = uitab(obj.htabgroup, 'Title', 'Result')
            obj.htab2 = uitab(obj.htabgroup, 'Title', 'Used Spots')
            obj.htabgroup2 = uitabgroup(obj.htab2);
            
            obj.handels.parentMenu = uimenu('Label','Alignment');
            obj.spotDetectorCam1  = guiSpotDetector(obj.htabgroup2,dataObject,4,nFrames,'Camera 1',obj.handels.parentMenu);
            obj.deletetry(obj.spotDetectorCam1.handels.parentMenu)
            obj.spotDetectorCam2  = guiSpotDetector(obj.htabgroup2,dataObject,3,nFrames,'Camera 2',obj.handels.parentMenu);
            obj.deletetry(obj.spotDetectorCam2.handels.parentMenu)
            
            uimenu(obj.handels.parentMenu,'Label','Estimate Alignment','Callback',@obj.alignmentEst);
            uimenu(obj.handels.parentMenu,'Label','Calculate Overlay','Separator','on','Callback',@obj.calculateOverlay);
            uimenu(obj.handels.parentMenu,'Label','Calculate Sub-Overlay','Callback',@obj.calculateSubOverlay);
            uimenu(obj.handels.parentMenu,'Label','Load','Separator','on','Callback',@obj.loadvars);
            uimenu(obj.handels.parentMenu,'Label','Save','Callback',@obj.savevars);
            
            obj.Figure = getParentFigure(h);  
            
            obj.paramsPreFilterFits = getDefaultParamsPreFilterFits;
            obj.paramsPreFilterFits.minPixelDist=2;
            obj.paramsPreFilterFits.clusterSizeMax=50;
            obj.paramsPreFilterFits.circularityMax=2;
            obj.spotDetectorCam1.paramsPreFilterFits = obj.paramsPreFilterFits;
            obj.spotDetectorCam2.paramsPreFilterFits = obj.paramsPreFilterFits;
            obj.paramsFit = getDefaultParamsFit;
            obj.paramsFit.FitSigma=true;

            obj.spotDetectorCam1.paramsFit=obj.paramsFit;
            obj.spotDetectorCam2.paramsFit=obj.paramsFit;

            obj.paramsFilterFits = getDefaultParamsFilterFits;
            obj.paramsFilterFits.minPixelDist=2;
            obj.spotDetectorCam1.paramsFilterFits = obj.paramsFilterFits;
            obj.spotDetectorCam2.paramsFilterFits = obj.paramsFilterFits;
      end
      
      function img = calculateOverlay(obj,C,~)
            if ~isempty(squeeze(obj.disttFormEst))
                
                if size(C,1) == 2 && size(C,2) == 2
                    pixelx = C(1,2):C(1,2)+C(2,2);
                    pixely = C(1,1):C(1,1)+C(2,1);
                else                
                    pixelx = [];
                    pixely=[];
                end
                gimr = zeros(size(obj.dataObject.getData(1,pixelx,pixely,[])));
                gimg = zeros(size(gimr));
                for i=1:size(gimg,3)
                    [image_out] = double(affine_trans( obj.dataObject.getData(1,pixelx,pixely,i), obj.initialtransVec(1:2), obj.initialtransVec(3:4), obj.initialtransVec(5)));
                    gimr(:,:,i) = imadjust(mat2gray(double(stretch(squeeze(obj.dataObject.getData(2,pixelx,pixely,i))))));                  
                    gimg(:,:,i) = imadjust(mat2gray(double(stretch(squeeze(image_out)))));             
                end
               
                gimg = permute(uint8(gimg * 256),[1 2 4 3]);
                color = obj.dataObject.getColor(1);
                imgg1 = cat(3,color(1).*gimg,color(2).*gimg,color(3).*gimg);
               
                gimr = permute(uint8(gimr * 256),[1 2 4 3]);
                color = obj.dataObject.getColor(2);                
                imgg2 = cat(3,color(1).*gimr,color(2).*gimr,color(3).*gimr);
                h = figure;
                img = imgg1+imgg2; 

                h2 = imtool3D(h,img);
            else 
                img = [];
            end
      end      
      function img = calculateSubOverlay(obj,C,~)
           
            if ~(size(C,1) == 2 && size(C,2) == 2)
                msize = size(obj.dataObject.getData(2),3);
                h = dipshow(obj.dataObject.getData(2,[],[],floor(linspace(1,msize,25))),'lin');
                [~,C] = dipcrop(h);
                close(h)
            end
            
            if size(C,1) == 2 && size(C,2) == 2  
                img = obj.calculateOverlay(C);   
            else
                img = [];
            end
            
      end  
      
        function dataChange(obj,src,~)
            notify(obj,'ParamChange');
        end
       
        function vs = savevars(obj,src,~)
            vs.class = obj.class;
            vs.tformTotal = obj.tformTotal;
            vs.rgbImage = obj.rgbImage;
            vs.disttFormEst = obj.disttFormEst;
            vs.cpPointsMoved = obj.cpPointsMoved;
            vs.cpPoints = obj.cpPoints;
            vs.idx = obj.idx;
            
            vs.paramsPreFilterFits = obj.paramsPreFilterFits;
            vs.paramsFit = obj.paramsFit;
            vs.paramsFilterFits = obj.paramsFilterFits;
            vs.spotDetectorCam1vs = obj.spotDetectorCam1.savevars;
            vs.spotDetectorCam2vs = obj.spotDetectorCam2.savevars;
            
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
                error(['This is no ' class ' class'])
            else
                obj.tformTotal = ws.tformTotal;
                obj.rgbImage = ws.rgbImage;
                obj.disttFormEst = ws.disttFormEst;
                obj.cpPointsMoved = ws.cpPointsMoved;
                obj.cpPoints = ws.cpPoints;
                obj.idx = ws.idx;
                obj.paramsPreFilterFits = ws.paramsPreFilterFits;
                obj.paramsFit = ws.paramsFit;
                obj.paramsFilterFits = ws.paramsFilterFits;
                obj.spotDetectorCam1.loadvars(ws.spotDetectorCam1vs,updateplot);
                obj.spotDetectorCam2.loadvars(ws.spotDetectorCam2vs,updateplot);
                obj.mseFormEst =  mean(obj.disttFormEst(obj.disttFormEst<1).^2)
                
                if updateplot
                    obj.plotImages;
                end
                obj.status=1;
            end   
        end
        
        
         function  savetoworkspace(obj,src,~)
                [file, pathname] = uigetfile([obj.class '.mat'] );
                if file ~= 0
                    answer{1} = fullfile(pathname, file);
                else
                    answer=[];
                end

            obj.whoAmI;
            evalin('base',[answer{1} ' = ' obj.name '.tformTotal']);           
         end
        
         function calculateImages(obj,src,~)
              if ~isempty(squeeze(obj.disttFormEst))
                    [image_out] = affine_trans(obj.mGellerCam1, obj.initialtransVec(1:2), obj.initialtransVec(3:4), obj.initialtransVec(5));
                    
                    gimr = imadjust(mat2gray(double(stretch(squeeze(obj.mGellerCam2))))); 
                    gimg = imadjust(mat2gray(double(stretch(squeeze(image_out))))); 
                    gimg = uint8(gimg * 256);
                    gimr = uint8(gimr * 256);

                    color = obj.dataObject.getColor(1);
                    imgg1 = cat(3,color(1).*gimg,color(2).*gimg,color(3).*gimg);

                    color = obj.dataObject.getColor(2);
                    imgg2 = cat(3,color(1).*gimr,color(2).*gimr,color(3).*gimr);
                    obj.rgbImage = imgg1+imgg2; 
              end
         end
       

        function plotImages(obj,src,~)
             if ~isempty(squeeze(obj.disttFormEst))
                    delete(obj.img1)
                    delete(obj.img2)
                    x1=0;
                    y1=0;

                    set(0,'CurrentFigure',obj.Figure);
                    sfh = subplot(1,2,1,'Parent', obj.htab1)
                    obj.img1 = scatter3(obj.cpPointsMoved(obj.disttFormEst<1,1),obj.cpPointsMoved(obj.disttFormEst<1,2),...
                    obj.cpPointsMoved(obj.disttFormEst<1,3),20,obj.disttFormEst(obj.disttFormEst<1),'filled')
                    colormap(sfh,jet);
                    caxis manual
                    caxis([min(obj.disttFormEst(obj.disttFormEst<1)) max(obj.disttFormEst(obj.disttFormEst<1))]);
                    hp4 = get(subplot(1,2,1),'Position')
                    colorbar('Position', [hp4(1)+hp4(3)  hp4(2)  0.01  hp4(2)+hp4(3)*2.1])
                    axis tight
                    shading interp;
                    view(90,90)
                    title(sprintf('Experimental mean square error: %0.2g\n',obj.mseFormEst))
           
                    
                    set(0,'CurrentFigure',obj.Figure);
                    subplot(1,2,2,'Parent', obj.htab1)
                    obj.img2 = subimage(obj.rgbImage);
                    xlim([0 size(obj.rgbImage,2)])
                    ylim([0 size(obj.rgbImage,1)])
                    hold on
                    
                    scatter(obj.cpPoints(obj.idx(obj.disttFormEst<1),1)+1,obj.cpPoints(obj.idx(obj.disttFormEst<1),2)+1,'xb')
                    cpPointsUn = transformPointsInverse(obj.tformTotal,obj.cpPointsMoved(:,1:2));
                    scatter(cpPointsUn(obj.disttFormEst<1,1)+1,cpPointsUn(obj.disttFormEst<1,2)+1,'or')                   
            end    
        end

        function  alignmentEst(obj,src,~)
            if obj.dataObject.status

                
                if isempty(obj.spotDetectorCam1.detParCam1) || isempty(obj.spotDetectorCam1.detParCam1)
                    obj.spotDetectorCam1.detectSpots;
                    obj.spotDetectorCam2.detectSpots;
                end
                
                obj.mGellerCam1 = mean(obj.dataObject.getData(3),3);
                obj.mGellerCam2 = mean(obj.dataObject.getData(4),3);
                
                [zoom, trans, ang] = fmmatch(obj.mGellerCam1,obj.mGellerCam2)
                if mean(trans) == 0 | mean(abs(trans)) > 10
                    sv2 = findshift(obj.mGellerCam1,obj.mGellerCam2,'iter');
                    [obj.initialtransVec,~] = find_affine_trans(obj.mGellerCam1, obj.mGellerCam2, [[1,1], sv2', 0]);
                else
                    [obj.initialtransVec,~] = find_affine_trans(obj.mGellerCam1, obj.mGellerCam2, [[1,1].*zoom, trans, ang]);
                end
                [~,R] = affine_trans(obj.mGellerCam1, obj.initialtransVec(1:2), obj.initialtransVec(3:4), obj.initialtransVec(5));

                maskCam1 = obj.spotDetectorCam1.maskFilt1&obj.spotDetectorCam1.maskPreFiltCam1;
                maskCam2 = obj.spotDetectorCam2.maskFilt1&obj.spotDetectorCam2.maskPreFiltCam1;
                
                obj.cpPoints = [obj.spotDetectorCam1.rawInitialFitResultsCam1.Coord(maskCam1,1) obj.spotDetectorCam1.rawInitialFitResultsCam1.Coord(maskCam1,2) obj.spotDetectorCam1.rawInitialFitResultsCam1.Frame(maskCam1)];
                obj.cpPointsMoved  = [obj.spotDetectorCam2.rawInitialFitResultsCam1.Coord(maskCam2,1) obj.spotDetectorCam2.rawInitialFitResultsCam1.Coord(maskCam2,2) obj.spotDetectorCam2.rawInitialFitResultsCam1.Frame(maskCam2)];
                obj.tformNoised = affine2d(R');
                  
                obj.cpPointsUnMoved = transformPointsInverse(obj.tformNoised,obj.cpPointsMoved(:,1:2));

                [obj.idx, obj.disttFormEst] = knnsearch(obj.cpPoints(:,1:2),obj.cpPointsUnMoved);
                a = pinv([obj.cpPoints(obj.idx(obj.disttFormEst<1),1:2) ones(size(obj.cpPoints(obj.idx(obj.disttFormEst<1),:) ,1),1)])* [obj.cpPointsUnMoved(obj.disttFormEst<1,:) ones(size(obj.cpPointsUnMoved(obj.disttFormEst<1,:),1),1)];

                Total = a*(R');
                Total(:,3)=[0;0;1;];
                obj.tformTotal = affine2d(Total);
                obj.cpPointsUnMoved = transformPointsInverse(obj.tformTotal,obj.cpPointsMoved(:,1:2));
                
                [obj.idx, obj.disttFormEst] = knnsearch(obj.cpPoints(:,1:2),obj.cpPointsUnMoved);              
                obj.mseFormEst =  mean(obj.disttFormEst(obj.disttFormEst<1).^2);
                
                obj.calculateImages;
                obj.plotImages;
                obj.status=1;
            else
                error('No bead data is loaded!')
            end
        end
        
        function rs = getSummary(obj)
              
            if isempty(obj.disttFormEst)
                rs = {'Estimate alignment first'};
            else
                 txt1 = sprintf('Experimental mean square error: %0.2g\n', mean(obj.disttFormEst.^2));
                rs = {txt1};
            end
        end
        
        function clear(obj,src,~)
            obj.paramsPreFilterFits=[];
            obj.paramsFit=[];
            obj.paramsFilterFits=[];
            obj.mGellerCam1=[];
            obj.mGellerCam2=[];
            obj.cpPoints=[];
            obj.tformNoised=[];
            obj.cpPointsMoved=[];
            obj.cpPointsUnMoved=[];
            obj.tformTotal=[];
            obj.disttFormEst=[];
            obj.idx=[];
            obj.status=0;
              
            delete(obj.img1)
            delete(obj.img2)

        end
        
   end
end