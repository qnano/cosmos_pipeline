classdef guiDriftCorr < handle
   properties
    %handles
    Figure
    Axes
    CLimit
    dataObject
    htab1
    imageStack
    
    %settings
    displayRGB
    cameraIndex
    
    %results
    sv
    drift
    C
    cutRegion
    RBGshift
    handels

    status=0;
    class = 'guiDriftCorr'
   end
   
    events
       ParamChange
   end
   
   methods
        function deletetry(obj,handle)
            if ~isempty(handle)
                try
                    delete(handle)
                catch
                end
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
      
        
      function obj = guiDriftCorr(h,dataObject,cameraIndex)
        if nargin < 3
            cameraIndex=1;
        end
        obj.status=0;
        
        obj.cameraIndex=cameraIndex;
        obj.dataObject = dataObject;
        obj.htab1 = uitab(h, 'Title', 'Drift');
        obj.handels.parentMenu = uimenu('Label','Drift');
        uimenu(obj.handels.parentMenu,'Label','Estimate Drift','Callback',@obj.driftEst);
        uimenu(obj.handels.parentMenu,'Label','Calculate Error RGB','Callback',@obj.calculateErrorRGB);
        uimenu(obj.handels.parentMenu,'Label','Load','Separator','on','Callback',@obj.loadvars);       
        uimenu(obj.handels.parentMenu,'Label','Save','Callback',@obj.savevars);            
      end               

       function vs = savevars(obj,src,~)
         
            vs.class = obj.class;
            vs.drift = obj.drift;
            vs.displayRGB = obj.displayRGB;
            vs.sv = obj.sv;
            vs.RBGshift = obj.RBGshift;
            vs.cutRegion = obj.cutRegion;
            vs.cameraIndex=obj.cameraIndex;
           
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
                
			
            class = obj.class;
            if ischar(answer{1})
                ws = load(answer{1});
                ws=ws.vs;
            elseif isstruct(src)
                ws=src;
            end
            if ~strcmpi(ws.class,class)
                error(['This is no ' class ' class'])
            else
                obj.drift=ws.drift;
                obj.displayRGB = ws.displayRGB;
                obj.sv = ws.sv;
                obj.RBGshift = ws.RBGshift;
                obj.cutRegion = ws.cutRegion;
                obj.cameraIndex=ws.cameraIndex;
                
                if updateplot
                    obj.plotImages;
                end
            end   
            obj.status=1;
        end

       function dataChange(obj,src,~)
           notify(obj,'ParamChange');
       end   

        function plotImages(obj,src,~)
            
            if ~isempty(squeeze(obj.sv))
                obj.deletetry( obj.imageStack)
                x1=0.05;
                y1=0.1;
                hAxis2 = subplot(1,2,1,'Parent', obj.htab1);
                if obj.displayRGB
                    obj.imageStack = imtool3D(hAxis2,obj.RBGshift);               
                else
                    color = obj.dataObject.getColor(obj.cameraIndex);
                    dat = permute(obj.cutRegion/255,[1 2 4 3]);
                    obj.imageStack = imtool3D(hAxis2,cat(3,color(1).*dat,color(2).*dat,color(3).*dat));                    
                end

                hAxis2 = subplot(1,2,2,'Parent', obj.htab1);
                plot(hAxis2,cumsum(obj.sv'),'linewidth',2)
                legend('x-drift','y-drift')

                ylabel('Cumultive Drift [pixel]')
                xlabel('t [frames]')
                grid on

                obj.drift=cumsum(obj.sv'); obj.drift=[[0 0]; obj.drift];
            end
        end

        function  driftEst(obj,C,~)
            if ~obj.dataObject.status
               error('No data is loaded for drift estimation!') 
            end
            
            %% plot pixel drift shift
            obj.displayRGB=false;
            if size(C,1) == 2 && size(C,2) == 2
                obj.C=C;
            else                
                msize = size(obj.dataObject.getData(obj.cameraIndex),3);
                color = obj.dataObject.getColor(obj.cameraIndex);
                data = obj.dataObject.getData(obj.cameraIndex,[],[],floor(linspace(1,msize,25)));
                h = dipshow(joinchannels('RGB',color(1).*data,color(2).*data,color(3).*data));
                [~,obj.C] = dipcrop(h);
                close(h)
            end
           
            obj.cutRegion = double(stretch( obj.dataObject.getData(obj.cameraIndex,obj.C(1,2):obj.C(1,2)+obj.C(2,2),obj.C(1,1):obj.C(1,1)+obj.C(2,1),[])));

            for i=1:size(obj.cutRegion,3)-1
                obj.sv(:,i) = findshift(obj.cutRegion(:,:,i),obj.cutRegion(:,:,i+1),'iter');
            end
            obj.sv=obj.sv;

            obj.plotImages;   
            obj.status=1;
        end
        
        function calculateErrorRGB(obj,src,~)
             obj.displayRGB=true;
             for i=1:size(obj.cutRegion,3)-1
                obj.sv(:,i) = findshift(obj.cutRegion(:,:,i),obj.cutRegion(:,:,i+1),'iter');
                if obj.displayRGB    
                    gimr = imadjust(mat2gray(double(stretch(squeeze( obj.cutRegion(:,:,i)))))); 
                    gimg = imadjust(mat2gray(double(stretch(squeeze(shift(obj.cutRegion(:,:,i+1),obj.sv(:,i))))))); 
                    gimg = uint8(gimg);
                    gimr = uint8(gimr);
                    filler = zeros(size(gimr),'uint8');
                    obj.RBGshift(:,:,:,i) = cat(3,gimr,gimg,filler);
                end
            end
            obj.sv=obj.sv;
            obj.plotImages;
        end
        
        function rs = getSummary(obj)             
            if isempty(obj.disttFormEst)
                rs = {'Estimate drift first'};
            else
                 txt1 = sprintf('Drift: mu = %0.2g std = %0.2g\n', mean(obj.drift),std(obj.drift));
                rs = {txt1};
            end
        end
        
        function clear(obj,src,~)
            obj.sv=[];
            obj.drift=[];
            
            obj.status=0;
        end
   end
end