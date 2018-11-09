classdef guiAlignment < handle
   properties
    Figure
    Axes
    CLimit
    dataObject
    htab1

    paramsPreFilterFits
    paramsFit
    paramsFilterFits      
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
        function  settingspff(obj,src,~)
            prompt={'maximum circularity','minimum  circularity','maximum PH1','minimum PH1','pixel distance','maximum clustersize','minimum clustersize'};
            name = 'params Pre-FilterFits';
            defaultans = {num2str(obj.paramsPreFilterFits.circularityMax), num2str(obj.paramsPreFilterFits.circularityMin),...
                num2str(obj.paramsPreFilterFits.PH1Max), num2str(obj.paramsPreFilterFits.PH1Min),...
                 num2str(obj.paramsPreFilterFits.minPixelDist), num2str(obj.paramsPreFilterFits.clusterSizeMin),...
                  num2str(obj.paramsPreFilterFits.clusterSizeMax)...
                  };
            options.Interpreter = 'tex';
            
            [answer] =  inputdlg(prompt,name,[1 50],defaultans,options);
            cancel = isempty(answer);
            if ~cancel
                obj.paramsPreFilterFits.circularityMax=str2num(answer{1});
                obj.paramsPreFilterFits.PH1Max=str2num(answer{2});
                obj.paramsPreFilterFits.PH1Min= str2num(answer{3});
                obj.paramsPreFilterFits.minPixelDist = str2num(answer{4});
                obj.paramsPreFilterFits.clusterSizeMin = str2num(answer{5});
                obj.paramsPreFilterFits.clusterSizeMax = str2num(answer{6});
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
        
      function obj = guiAlignment(h,dataObject)
            obj.dataObject = dataObject;
            obj.htab1 = uitab(h, 'Title', 'Alignment');
            f = uimenu('Label','Alignment');
            uimenu(f,'Label','Pre-FilterFit Settings','Callback',@obj.settingspff);
            uimenu(f,'Label','Fit Settings','Callback',@obj.settingsf);
            uimenu(f,'Label','FilterFit Settings','Callback',@obj.settingsff);
            uimenu(f,'Label','Estimate Alignment','Callback',@obj.alignmentEst);
            uimenu(f,'Label','Save','Callback',@obj.savevars);

            obj.paramsPreFilterFits = getDefaultParamsPreFilterFits;
            obj.paramsPreFilterFits.minPixelDist=2;
            obj.paramsPreFilterFits.clusterSizeMax=25;
            obj.paramsPreFilterFits.circularityMax=2;

            obj.paramsFit = getDefaultParamsFit;
            obj.paramsFit.FitSigma=true;

            obj.paramsFilterFits = getDefaultParamsFilterFits;
            obj.paramsFilterFits.minPixelDist=2;
            
            obj.Figure = getParentFigure(h);

            
      end
  
       function dataChange(obj,src,~)
           notify(obj,'ParamChange');
       end
        
        
         function  savevars(obj,src,~)
            prompt={'Variable name'};
            name = 'Save';
            defaultans = {'tformTotal'};
            options.Interpreter = 'tex';
            answer = inputdlg(prompt,name,[1 40],defaultans,options);

            obj.whoAmI;
            evalin('base',[answer{1} ' = ' obj.name '.tformTotal']);           
        end
       

        function plotImages(obj,src,~)

            if ~isempty(squeeze(obj.disttFormEst))
                    delete(obj.img1)
                    delete(obj.img2)
                    x1=0;
                    y1=0;
                    
                    set(0,'CurrentFigure',obj.Figure);
                    sfh = subplot(1,2,1,'Parent', obj.htab1)
                    obj.img1 = scatter3(obj.cpPointsMoved(obj.idx(obj.disttFormEst<1),1),obj.cpPointsMoved(obj.idx(obj.disttFormEst<1),2),...
                    obj.cpPointsMoved(obj.idx(obj.disttFormEst<1),3),20,obj.disttFormEst(obj.disttFormEst<1),'filled')
                    colormap(sfh,jet);
                    caxis manual
                    caxis([0 max(1,max(obj.disttFormEst))]);
                    colorbar;
                    shading interp;
                    title(sprintf('Experimental mean square error: %0.2g\n', mean(obj.disttFormEst.^2)))

                    [image_out] = affine_trans(obj.mGellerCam1, obj.initialtransVec(1:2), obj.initialtransVec(3:4), obj.initialtransVec(5));
                                      
                    gimr = mat2gray(double(stretch(squeeze(obj.mGellerCam2)))); 
                    gimg = mat2gray(double(stretch(squeeze(image_out)))); 
                    gimg = uint8(gimg * 256);
                    gimr = uint8(gimr * 256);
                    filler = zeros(size(gimr),'uint8');
                    obj.rgbImage = cat(3,gimr,gimg,filler);
                    
                    set(0,'CurrentFigure',obj.Figure);
                    subplot(1,2,2,'Parent', obj.htab1)
                    obj.img2 = subimage(obj.rgbImage);
                    
            end    
        end

        function  alignmentEst(obj,src,~)
            if ~isempty(obj.dataObject.beadDataCam1)

                [coordsCam1,detParCam1,~] = LLRMapv2(obj.dataObject.beadDataCam1,obj.paramsFit.PSFSigma);
                [coordsCam2,detParCam2,~] = LLRMapv2(obj.dataObject.beadDataCam2,obj.paramsFit.PSFSigma);

                [ maskPreFiltCam1 ] =  preFilterFits(coordsCam1,detParCam1,obj.paramsPreFilterFits);

                %% Pre Filter Detection Clusters Cam2

                [ maskPreFiltCam2 ] =  preFilterFits(coordsCam2,detParCam2,obj.paramsPreFilterFits);

                coodsUnCut=round(coordsCam1+(1.5*(2*obj.paramsFit.PSFSigma+1)-0.5).*[ones(size(coordsCam1,1),1) ones(size(coordsCam1,1),1) zeros(size(coordsCam1,1),1)]);
                [ rawFitResultsCam1 ] = fitBoxCenters( single(squeeze(obj.dataObject.beadDataCam1)),coodsUnCut,obj.paramsFit);
                [ maskFilt1 ] =  filterFits(rawFitResultsCam1,obj.paramsFilterFits);

                coodsUnCut=round(coordsCam2+(1.5*(2*obj.paramsFit.PSFSigma+1)-0.5).*[ones(size(coordsCam2,1),1) ones(size(coordsCam2,1),1) zeros(size(coordsCam2,1),1)]);
                [ rawFitResultsCam2 ] = fitBoxCenters( single(squeeze(obj.dataObject.beadDataCam2)),coodsUnCut,obj.paramsFit);

                [ maskFilt2 ] =  filterFits(rawFitResultsCam2,obj.paramsFilterFits);

                %%
                
                obj.mGellerCam1 = mean(obj.dataObject.beadDataCam1,3);
                obj.mGellerCam2 = mean(obj.dataObject.beadDataCam2,3);
                sv2 = findshift(obj.mGellerCam1,obj.mGellerCam2,'iter');

                [obj.initialtransVec,~] = find_affine_trans(obj.mGellerCam1, obj.mGellerCam2, [1, 1, sv2', 0]);

                [~,R] = affine_trans(obj.mGellerCam1, obj.initialtransVec(1:2), obj.initialtransVec(3:4), obj.initialtransVec(5));


                obj.cpPoints = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1) rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,2) rawFitResultsCam1.Frame(maskPreFiltCam1&maskFilt1)];
                obj.cpPointsMoved  = [rawFitResultsCam2.Coord(maskPreFiltCam2&maskFilt2,1) rawFitResultsCam2.Coord(maskPreFiltCam2&maskFilt2,2) rawFitResultsCam2.Frame(maskPreFiltCam2&maskFilt2)];
                obj.tformNoised = affine2d(R');
                obj.cpPointsUnMoved = transformPointsInverse(obj.tformNoised,obj.cpPoints(:,1:2));

                [obj.idx, obj.disttFormEst] = knnsearch(obj.cpPointsMoved(:,1:2),obj.cpPointsUnMoved);

                a = pinv([obj.cpPointsMoved(obj.idx(obj.disttFormEst<1),1:2) ones(size(obj.cpPointsMoved(obj.idx(obj.disttFormEst<1),:) ,1),1)])* [obj.cpPointsUnMoved(obj.disttFormEst<1,:) ones(size(obj.cpPointsUnMoved(obj.disttFormEst<1,:),1),1)];

                Total = a*(R');
                Total(:,3)=[0;0;1;];
                obj.tformTotal = affine2d(Total);
                obj.cpPointsUnMoved = transformPointsInverse(obj.tformTotal,obj.cpPoints(:,1:2));

                [~, obj.disttFormEst] = knnsearch(obj.cpPointsMoved(obj.idx(obj.disttFormEst<1),1:2),transformPointsInverse(obj.tformTotal,obj.cpPoints(obj.disttFormEst<1,1:2) ));
                
                obj.plotImages;
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
              
            delete(obj.img1)
            delete(obj.img2)

        end
        
   end
end