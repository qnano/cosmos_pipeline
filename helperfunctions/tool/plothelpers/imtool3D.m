classdef imtool3D < handle
    properties
    parent
    rgbImage
    himg
    haxis
    parentFigure
    Slider
    stackLength
    pStretch
    l
    hdragline
    
    showSliceMethod
    params
    imgIntensityMapping
    noSlider
    
    eventCount
    end
    
    events
        filterChange
    end
    
    methods
        
        function obj = imtool3D(parent,imgrgb,params)
            if nargin < 3
                params=[];
            end
            obj.params=params;
            obj.pStretch.low = 0;    
            obj.pStretch.high = 1;
            obj.noSlider = false;
            obj.eventCount=0;
%             obj.pStretch.minimum = 0;
%             obj.pStretch.maximum = 1;
            
            if isa(imgrgb,'function_handle')
                obj.showSliceMethod=imgrgb;
            else
                obj.rgbImage=imgrgb;
                if size(imgrgb,3) == 1
                    obj.noSlider=true;
                end
            end
            
            obj.imgIntensityMapping='linear';
            
            if isa(parent,'matlab.graphics.axis.Axes')
                obj.parent=get(parent,'Parent');
                obj.haxis=parent;                
            else
                obj.parent=parent;
                obj.haxis= subplot(1,1,1,'Parent', obj.parent);
            end          
            
            if ~isempty(obj.showSliceMethod)
                [img, obj.stackLength] = obj.showSliceMethod(obj.params);
                obj.himg= imshow(img,'parent',obj.haxis);
                obj.showSliceMethod(obj,obj.himg,1);
            else
                if ndims(obj.rgbImage) <= 3
                    a = sort(reshape(obj.rgbImage(:,:,1),[1 size(obj.rgbImage,1)*size(obj.rgbImage,2)]));
                    img = imadjust(mat2gray(obj.rgbImage(:,:,1)),[  a(max(1,ceil(size(a,2)*obj.pStretch.low)))./max(a) a(min(end,round(size(a,2)*obj.pStretch.high)))./max(a)],[0 1]);
                elseif ndims(obj.rgbImage) == 4 && size(obj.rgbImage,3) == 3
                    img = obj.rgbImage(:,:,:,1);
                end 
                obj.himg= imshow(img,'parent',obj.haxis);
                temp = size(obj.rgbImage);
                obj.stackLength = temp(end);
                title(sprintf('t=%d',round(1)),'parent',obj.haxis)
            end
            axis off
            
            obj.parentFigure = obj.getParentFigure(obj.parent); 
            
            c = uicontextmenu('Parent', obj.parentFigure );
            obj.himg.UIContextMenu=c;
             
            uimenu(c,'Label','linear','Callback',@obj.setlinestyle);
            uimenu(c,'Label','log','Callback',@obj.setlinestyle);  
            uimenu(c,'Label','limits','Callback',@obj.setlinestyle);
            if ~isa(obj.rgbImage,'function_handle')
                uimenu(c,'Label','save as tiff','Separator','on','Callback',@obj.saveAsTiff);
            end
            
            pos=obj.himg.Parent.Position;
            Newpos= [pos(1)-0.06 pos(2) 0.05 pos(4)];
            
            if ~obj.noSlider
                obj.Slider=uicontrol('style','slider',...
                  'units','normalized','position',Newpos,...
                  'min',1,'max',obj.stackLength,'value',1,'Parent', obj.parent); %,'Callback',@obj.sliderCallBack
                obj.Slider.SliderStep=[1 ceil(obj.stackLength*0.1)]./obj.Slider.Max;
                addlistener(obj.Slider, 'ContinuousValueChange', @obj.sliderCallBack);
            end

  
            if ~isfield(obj.parentFigure.UserData,'scrollEventClass')
                obj.parentFigure.UserData.scrollEventClass = scrollEvent(obj.parentFigure);
            end
            addlistener(obj.parentFigure.UserData.scrollEventClass,'scroll',@obj.imgScroll);
            addlistener(obj,'filterChange',@obj.filterUpdate);
            
        end
        
        function saveAsTiff(obj,src,~)
                [file, pathname] = uiputfile(['data' date '.tiff'] );
                if file ~= 0
                    answer{1} = fullfile(pathname, file);
                else
                    answer=[];
                end
               
%                 metadata = createMinimalOMEXMLMetadata(uint8(obj.rgbImage));

                if ~isempty(answer)
                    bfsave(obj.rgbImage,answer{1},'dimensionOrder', 'XYCZT'); %,'BigTiff', true
                end
        end
        
        function filterUpdate(obj,source,~)
            obj.showSlice(round(get(obj.Slider,'value')));
        end
        
        function setlinestyle(obj,source,~)
            switch source.Label
                case 'linear'
                    obj.imgIntensityMapping='linear';
                case 'log'
                    obj.imgIntensityMapping='log';
                case 'limits'
                    figure; 
                    hsub1=subplot(1,1,1);
                    if ~isempty(obj.showSliceMethod)
                        img = obj.showSliceMethod(obj,obj.himg,round(get(obj.Slider,'value')),false);
                        obj.rgbImage=img(:,:,1);
                    end
                    
                    hist(hsub1,reshape(obj.rgbImage(:,:,1),[1 size(obj.rgbImage,1)*size(obj.rgbImage,2)]),size(obj.rgbImage,1))
                    obj.l(1)=line([1 1].*max(max(obj.rgbImage(:,:,1))),get(hsub1, 'YLim'),'Color','r','Linewidth',2,'parent',hsub1,'tag','line2','hittest','off');
                    obj.l(2)=line([1 1].*min(min(obj.rgbImage(:,:,1))),get(hsub1, 'YLim'),'Color','g','Linewidth',2,'parent',hsub1,'tag','line1','hittest','off'); %min(obj.rawInitialFitResultsCam1.Photons)
                    set(hsub1 ,'buttondownfcn',@obj.clickline,'xlimmode','manual');       
%                     prompt={'low procent','high procent'};
%                     name = 'Enter image limits';
%                     defaultans = {num2str(obj.pStretch.low),num2str(obj.pStretch.high)};
%                     options.Interpreter = 'tex';
%                     answer = inputdlg(prompt,name,[1 40],defaultans,options);
%                     obj.pStretch.low=str2double(answer{1});
%                     obj.pStretch.high=str2double(answer{2});
                    %obj.pStretch.minimum=str2double(answer{3});
                    %obj.pStretch.maximum=str2double(answer{4});
            end
       
            obj.showSlice(round(get(obj.Slider,'value')));
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
            
            xdata = get(obj.l(1),'xdata');
            obj.pStretch.high=max(0,xdata(1)./max(max(obj.rgbImage(:,:,1))));
            xdata = get(obj.l(2),'xdata');
            obj.pStretch.low=min(1,xdata(1)./max(max(obj.rgbImage(:,:,1))));
            notify(obj,'filterChange');
        end
        
        function sliderCallBack(obj,~,~)
            newSice=round(get(obj.Slider,'value'));
            obj.showSlice(newSice);
        end
        
        
        function imgScroll(obj,~,~,~)
            if obj.isMouseOverAxes(obj.haxis)
                newSlice=round(get(obj.Slider,'value')-obj.parentFigure.UserData.scrollEventClass.VerticalScrollCount);
                if ~isempty(newSlice) && newSlice>=1 && obj.stackLength >= newSlice
                    set(obj.Slider,'value',newSlice);
%                     obj.eventCount=obj.eventCount+1;
%                     notify(obj,'filterChange');
%                     obj.eventCount=obj.eventCount-1;
%                     obj.showSlice(newSlice);     
%                 notify(obj.Slider,'ContinuousValueChange')
                end
            end
        end
        
        function showSlice(obj,newSlice)
            obj.eventCount=obj.eventCount+1;
            if isempty(newSlice)
                newSlice=1;
            end         
            if obj.eventCount ==1
                if ~isempty(obj.showSliceMethod)
                    obj.showSliceMethod(obj,obj.himg,newSlice);
                else

                    if ndims(obj.rgbImage) <= 3
                        switch (obj.imgIntensityMapping)
                            case 'linear'
                                a = sort(reshape(obj.rgbImage(:,:,newSlice),[1 size(obj.rgbImage,1)*size(obj.rgbImage,2)]));
                                obj.himg.CData = double(imadjust(mat2gray((obj.rgbImage(:,:,newSlice))),[  a(max(1,ceil(size(a,2)*obj.pStretch.low)))./max(a) a(min(end,round(size(a,2)*obj.pStretch.high)))./max(a)],[0 1]));
                            case 'log'
                                imglog = (real(log(1+obj.rgbImage(:,:,newSlice))));
    %                             a = sort(reshape(imglog,[1 size(obj.rgbImage,1)*size(obj.rgbImage,1)]));
                                obj.himg.CData = mat2gray(double(stretch(imglog,[  obj.pStretch.low],[ obj.pStretch.high*100],[1],[256])));
                        end
                    elseif ndims(obj.rgbImage) == 4 && size(obj.rgbImage,3) == 3
                         switch (obj.imgIntensityMapping)
                            case 'linear'
                                a = sort(reshape(obj.rgbImage(:,:,1,newSlice),[1 size(obj.rgbImage,1)*size(obj.rgbImage,2)]));
    %                             a2 = sort(reshape(obj.rgbImage(:,:,2,newSlice),[1 size(obj.rgbImage,1)*size(obj.rgbImage,2)]));
    %                             a3 = sort(reshape(obj.rgbImage(:,:,3,newSlice),[1 size(obj.rgbImage,1)*size(obj.rgbImage,2)]));
                                obj.himg.CData = cat(3,...
                                    double(imadjust(mat2gray((obj.rgbImage(:,:,1,newSlice))),double([  a(max(1,ceil(size(a,2)*obj.pStretch.low)))./max(a) a(min(end,round(size(a,2)*obj.pStretch.high)))./max(a)]),round([0 1]))),...
                                    double(imadjust(mat2gray((obj.rgbImage(:,:,2,newSlice))),double([  a(max(1,ceil(size(a,2)*obj.pStretch.low)))./max(a) a(min(end,round(size(a,2)*obj.pStretch.high)))./max(a)]),[0 1])),...
                                    double(imadjust(mat2gray((obj.rgbImage(:,:,3,newSlice))),double([  a(max(1,ceil(size(a,2)*obj.pStretch.low)))./max(a) a(min(end,round(size(a,2)*obj.pStretch.high)))./max(a)]),[0 1])));
                            case 'log'
                                imglog1 = (real(log(1+obj.rgbImage(:,:,1,newSlice))));
                                imglog1 = mat2gray(double((obj.rgbImage(:,:,1,newSlice)>0).*stretch(imglog1,[  obj.pStretch.low],[ obj.pStretch.high*100],[1],[256])));
                                imglog2 = (real(log(1+obj.rgbImage(:,:,2,newSlice))));
                                imglog2 = mat2gray(double((obj.rgbImage(:,:,2,newSlice)>0).*stretch(imglog2,[  obj.pStretch.low],[ obj.pStretch.high*100],[1],[256])));
                                imglog3 = (real(log(1+obj.rgbImage(:,:,3,newSlice))));
                                imglog3 = mat2gray(double((obj.rgbImage(:,:,3,newSlice)>0).*stretch(imglog3,[  obj.pStretch.low],[ obj.pStretch.high*100],[1],[256])));
        %                     a = sort(reshape(imglog,[1 size(obj.rgbImage,1)*size(obj.rgbImage,1)]));
                                obj.himg.CData = cat(3,imglog1,imglog2,imglog3);
        %                     obj.himg.CData = obj.rgbImage(:,:,:,newSlice);
                        end
                    end 
                    title(sprintf('t=%d',round(newSlice)),'parent',obj.haxis)
                end
            else
                disp('nope')
            end
            obj.eventCount=obj.eventCount-1;
        end
        
        function overAxes = isMouseOverAxes(~,ha)
            %This function checks if the mouse is currently hovering over the axis in
            %question. hf is the handle to the figure, ha is the handle to the axes.
            %This code allows the axes to be embedded in any size heirarchy of
            %uipanels.
           
            point = get(ha,'CurrentPoint');
            x = point(1,1); y = point(1,2);
            xlims = get(ha,'Xlim');
            ylims = get(ha,'Ylim');

            overAxes = x>=xlims(1) & x<=xlims(2) & y>=ylims(1) & y<=ylims(2);
        end
        
        function fig = getParentFigure(~,fig)
            while ~isempty(fig) && ~strcmp('figure', get(fig,'type'))
              fig = get(fig,'parent');
            end
        end
    end
   
end




