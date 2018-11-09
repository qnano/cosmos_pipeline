classdef guiHMMAnalyse < guiAnalyzeDynamics
    properties
    initialClass
    segModel
    sigModel
    dataTransSignal
    indexInital
    evidence
    prior
    NTimes = 2;
    NOrder=5;
    fittedTraces;
    selectedModel
    mdSliderValue
    mgSliderValue
    maximumlikelihood = 0;
    
    sortState=1;
    end
    methods
        
        function obj = guiHMMAnalyse(h,dataObject,alignmentObject,driftObject,spotDetectorObject,intergrationTime,menuhandle,menuName,startFrame,detectionFrame,colocalized,flip)
             if nargin < 8 || isempty(menuName)
                menuName = '';
            end
            
            if nargin <  6 || isempty(intergrationTime )
                intergrationTime = 1;
            end
            if nargin < 7 || isempty(menuhandle)
                menuhandle=[];
            end
            if nargin < 8 || isempty(menuName)
            	menuName='HMM';
            else
                menuName= ['HMM' menuName];
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

             if isa(spotDetectorObject,'spotDetectorMultiplexRNAWrapper')
                spotDetectorObject = spotDetectorObject.spotDetectorObjects{end};
             end
             
            obj@guiAnalyzeDynamics(h,dataObject,alignmentObject,driftObject,spotDetectorObject,intergrationTime,menuhandle,menuName,startFrame,detectionFrame,colocalized,flip);
            
            obj.classSingle = 'guiHMMw';
            obj.classAdded = 'guiHMMAdded';
            
            obj.deletetry(obj.handels.childMenu);
            
            obj.handels.childMenu(1) = uimenu('Parent',obj.handels.parentMenu,'Label','Select Subregion','Callback',@obj.selectSubRegion);
            obj.handels.childMenu(2) = uimenu('Parent',obj.handels.parentMenu,'Label','Threshold Settings','Callback',@obj.setParamsGeneral);
            obj.handels.childMenu(3) = uimenu('Parent',obj.handels.parentMenu,'Label','Rastergram Settings','Callback',@obj.setRastergramSettings)
            obj.handels.childCheckMenu(1) = uimenu('Parent',obj.handels.parentMenu,'Label','Analyze','Separator','on','Callback',@obj.analyze);
            obj.handels.childCheckMenu(2) = uimenu('Parent',obj.handels.parentMenu,'Label','Variational Bayes HMM Signal Fit','Callback',@obj.bayesSignalFitGMM);

            obj.handels.childMenu(10) = uimenu('Parent',obj.handels.parentMenu,'Label','Animation','Callback',@obj.exportAnimation);
            obj.handels.childMenu(5) = uimenu('Parent',obj.handels.parentMenu,'Label','Calculate Rastergram','Callback',@obj.runRastergram);
            obj.handels.childMenu(6) = uimenu('Parent',obj.handels.parentMenu,'Label','Bootstrap Rastergram','Callback',@obj.bootStrapRastergramMenu);
            obj.handels.childMenu(11) = uimenu('Parent',obj.handels.parentMenu,'Label','Export Data','Callback',@obj.exportData);
            obj.handels.childMenu(7) = uimenu('Parent',obj.handels.parentMenu,'Label','Load','Separator','on','Callback',@obj.loadvars);  
            obj.handels.childMenu(8) = uimenu('Parent',obj.handels.parentMenu,'Label','Add','Callback',@obj.addvars);              
            obj.handels.childMenu(9) = uimenu('Parent',obj.handels.parentMenu,'Label','Save','Callback',@obj.savevars);
            htabgroup = obj.handels.rastergramTap.Parent;
            
            slmin = 0;
            slmax = 50;
             obj.handels.mdSlider = uicontrol('Parent',  obj.handels.traceEstimateTab,'Style','Slider','SliderStep',[1 1]./(slmax-slmin)/2,...
            'Units','normalized','Position',[0.8 0.97 0.1 0.02],'Min',slmin,'Max',slmax,...
                  'Value',1); 
              
            obj.handels.mgSlider = uicontrol('Parent',  obj.handels.traceEstimateTab,'Style','Slider','SliderStep',[1 1]./(slmax-slmin)/2,...
            'Units','normalized','Position',[0.8 0.95 0.1 0.02],'Min',slmin,'Max',slmax,...
                  'Value',1); 
              
            set(obj.handels.mdSlider,'Callback',@obj.updateFilter);
            set(obj.handels.mgSlider,'Callback',@obj.updateFilter);
            
            obj.handels.Figure = getParentFigure(h);
            if ~isfield(obj.handels.Figure.UserData,'sliderFilterEventClass')
                obj.handels.Figure.UserData.sliderFilterEventClass = sliderEventClass(...
                    [get(obj.handels.mdSlider,'Value'); get(obj.handels.mgSlider,'Value')],...
                    [ get(obj.handels.mdSlider,'Max'); get(obj.handels.mgSlider,'Max')]);
            else
                TH = obj.handels.Figure.UserData.sliderFilterEventClass.TH;
                obj.handels.mgSlider.Value = TH(2);
                obj.handels.mdSlider.Value = TH(1);    
            end
            obj.mgSliderValue=obj.handels.mgSlider.Value;
            obj.mdSliderValue=obj.handels.mdSlider.Value;
            addlistener(obj.handels.Figure.UserData.sliderFilterEventClass,'sliderEvent',@obj.changeFilterPosition)
            
            
            obj.prior.alpha = 1;
            obj.prior.kappa = 100;
            obj.prior.m = '1';
            d=2; 
            obj.prior.v = 10;
            obj.prior.Ws = 10;
            obj.prior.alphaA=1;
            
            obj.handels.checkedMenuLine=[];
            obj.handels.fittedLine=[];
            
            obj.handels.filterSlider = uicontrol('Parent',  obj.handels.traceEstimateTab,'Style','text','Units','normalized','position',[0.9 0.91 0.1 0.08],'String',sprintf('Max gab = %d \n Min duration = %d', max(1,obj.handels.mgSlider.Value), max(1,obj.handels.mdSlider.Value)));

        end
        
        function runRastergram(obj,src,~)
             if ~obj.status
                error('Analysis has to be completed before rastergram can be calculated!')
             end
            if sum(obj.spotsIncluded) > 2
                [traces2_allfr, ttb] = obj.getTraceMatrix;

                obj.calculateRastergram(traces2_allfr,ttb,[],obj.sortState);
            else
                warning('Not enough traces to create rastergram.')
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

                    [meanDwellTimeBS,meanTime2FirstEventBS] = obj.bootStrapRastergram(obj.Nbootstraps,obj.sortState);
                end
              else
                  h = msgbox('Perform analysis first')
              end
        end
        
        function selectMenu(obj,index)
            for i=1:size(obj.handels.childCheckMenu,2)
                obj.handels.childCheckMenu(i).Checked = 'off';
            end
            obj.handels.childCheckMenu(index).Checked = 'on';
        end
        
        function idx = getSelectedMenu(obj)
            idx = zeros(size(obj.handels.childCheckMenu,2));
            for i=1:size(obj.handels.childCheckMenu,2)
                if strcmp(obj.handels.childCheckMenu(i).Checked,'on')
                    idx(i)=1;
                end
            end
        end
        
        function addLineHandels(obj,lineHandle)
            obj.handels.checkedMenuLine(end+1) = lineHandle;
        end
        
        function setLineHandels(obj,index,lineHandle)
            obj.handels.checkedMenuLine(index) = lineHandle;
        end
        
        function deleteLineHandels(obj)
            for i=1:size(obj.handels.checkedMenuLine)
                deletetry(obj.handels.checkedMenuLine(i))
            end
        end
        
        function updateFilter(obj,src,~)
            obj.handels.Figure.UserData.sliderFilterEventClass.changePosition( [get(obj.handels.mdSlider,'Value'); get(obj.handels.mgSlider,'Value')]);
        end
        

        function changeFilterPosition(obj,src,~)  
            obj.handels.mdSlider.Value = obj.handels.Figure.UserData.sliderFilterEventClass.TH(1);
            obj.handels.mgSlider.Value= obj.handels.Figure.UserData.sliderFilterEventClass.TH(2);
            obj.mdSliderValue = obj.handels.mdSlider.Value;
            obj.mgSliderValue = obj.handels.mgSlider.Value;
                obj.deleteLineHandels;
                if isfield(obj.handels,'selectedSpotIndex')
                    idx = find(obj.getSelectedMenu,1,'first');
                    if isempty(idx)
                        idx=1;
                        obj.selectMenu(1);
                    end
                    switch(idx)
                        case {1}
                            [traces2_allfr] = obj.getBinary(obj.handels.selectedSpotIndex);
                        case{2}
                            traces2_allfr = obj.fittedTraces(obj.handels.selectedSpotIndex,:);
                    end

                    se1 = strel('line', max(1,obj.handels.mdSlider.Value), 0);
                    se2 = strel('line', max(1,obj.handels.mgSlider.Value), 0);

                    A1 = imclose(traces2_allfr,se1);
                    initial = imopen(A1,se2);

                    ymax = max(max(obj.bg1I(obj.handels.selectedSpotIndex,:,2-0)*(obj.PSFSigma*2)^2),max(max(obj.photons1I(obj.handels.selectedSpotIndex,:,2-0)),max(obj.photons1I(obj.handels.selectedSpotIndex,:,1+0)))); 

                    if isfield(obj.handels,'filterdSegmentation') 
                        obj.deletetry(obj.handels.filterdSegmentation)
                    end

                    obj.handels.filterdSegmentation = line(1:size(obj.photons1I,2),double(initial)./max(max(double(initial))).*ymax...
                    ,'Color','k','parent',obj.handels.traceEstimateTabPlot,'tag','background','hittest','off');
                obj.handels.filterSlider.String = sprintf('Max gab = %d \n Min Duration = %d',round(max(1,obj.handels.mdSlider.Value)),round(max(1,obj.handels.mgSlider.Value)));
                end
        end        
        
        function plotBayesLines(obj,src,~)
            idx = find(obj.getSelectedMenu,1,'first');
            obj.deletetry(obj.handels.fittedLine)
            switch(idx)
                case {2,3}
                    if ~isempty(obj.fittedTraces)
                            ymax = max(max(obj.bg1I(obj.handels.selectedSpotIndex,:,2-0)*(obj.PSFSigma*2)^2),...
                            max(max(obj.photons1I(obj.handels.selectedSpotIndex,:,2-0)),...
                            max(obj.photons1I(obj.handels.selectedSpotIndex,:,1+0)))); 
                        obj.handels.fittedLine = line(1:size(obj.photons1I,2),obj.fittedTraces(obj.handels.selectedSpotIndex,:)./max(max(obj.fittedTraces)).*ymax...
                        ,'Color','k','parent',obj.handels.traceEstimateTabPlot,'tag','background','hittest','off');
                    end
                case {4,5}
                case 1
                otherwise
                    error('unknown line')
            end
        end
        
       function bayesSignalFitGMM(obj,src,~)
           
              if nargin > 1 && ischar(src)
                switch upper(src)
                    case 'MAXIMUMLIKELIHOOD'
                        obj.maximumlikelihood=1;
                    case 'POSTERIOR'
                        obj.maximumlikelihood=0;
                end
              end
        
            idx = find(obj.getSelectedMenu,1,'first');
%             if isempty(obj.fittedTraces) || idx == 1
                [traces2_allfr,data] = obj.getBinary(obj.spotsIncluded);

                data=data(:,:,3:4);
                se1 = strel('line',  max(1,obj.handels.mdSlider.Value), 0);
                se2 = strel('line',  max(1,obj.handels.mgSlider.Value), 0);

                A1 = imclose(traces2_allfr,se1);
                obj.initialClass = imopen(A1,se2);
                indexInital = sum(obj.initialClass,2)>1;

                idxData = repmat(obj.initialClass,[1 1 size(data,3)]);
                storeData = nan(size(data));
                storeData(idxData) = data(idxData);

                storeDataSelected = storeData(indexInital,:,:);
                dataTrans = permute(storeDataSelected,[3 2 1]);
                Title = 'VBHMM Settings';

                %%%% SETTING DIALOG OPTIONS
                Options.Resize = 'on';
                Options.Interpreter = 'tex';
                Options.CancelButton = 'on';
                Options.ApplyButton = 'off';

                Option.Dim = 1; % Horizontal dimension in fields

                Prompt = {};
                Formats = {};
                
                
                Prompt(1,:) = {'Maximum Order','NOrder',[]};
                Formats(1,1).type = 'edit';
                Formats(1,1).format = 'text';
                Formats(1,1).size = [-1 0];
                Formats(1,1).span = [1 1]; 
                DefAns.NOrder = mat2str(obj.NOrder);

                Prompt(2,:) = {'N Samples','NTimes',[]};
                Formats(2,1).type = 'edit';
                Formats(2,1).format = 'text';
                Formats(2,1).size = [-1 0];
                Formats(2,1).span = [1 1]; 
                DefAns.NTimes = mat2str(obj.NTimes);
                
                Prompt(3,:) = {'Prior \alpha','alpha',[]};
                Formats(3,1).type = 'edit';
                Formats(3,1).format = 'text';
                Formats(3,1).size = [-1 0];
                Formats(3,1).span = [1 1]; 
                DefAns.alpha = mat2str(obj.prior.alpha);

                Prompt(4,:) = {'Prior \kappa','kappa',[]};
                Formats(4,1).type = 'edit';
                Formats(4,1).format = 'text';
                Formats(4,1).size = [-1 0];
                Formats(4,1).span = [1 1]; 
                DefAns.kappa = mat2str(obj.prior.kappa);                

                Prompt(5,:) = {'Prior m','m',[]};
                Formats(5,1).type = 'edit';
                Formats(5,1).format = 'text';
                Formats(5,1).size = [-1 0];
                Formats(5,1).span = [1 1]; 
                DefAns.m = mat2str(obj.prior.m);                

                Prompt(6,:) = {'Prior v','v',[]};
                Formats(6,1).type = 'edit';
                Formats(6,1).format = 'text';
                Formats(6,1).size = [-1 0];
                Formats(6,1).span = [1 1]; 
                DefAns.v = mat2str(obj.prior.v);  

                Prompt(7,:) = {'Prior W^{1/2}','Ws',[]};
                Formats(7,1).type = 'edit';
                Formats(7,1).format = 'text';
                Formats(7,1).size = [-1 0];
                Formats(7,1).span = [1 1]; 
                DefAns.Ws = mat2str(obj.prior.Ws);                  

                Prompt(8,:) = {'Fit type','List',[]};
                Formats(8,1).type = 'list';
                Formats(8,1).style = 'popupmenu';
                Formats(8,1).items = {'Maximum likelihood','Posterior probability'};
                DefAns.List =obj.maximumlikelihood+1;
                
%                 prompt={'Maximum Order','','Prior \alpha','','Prior m','Prior v',''};
%                 name = 'Threshold Params';
% 
%                 defaultans = {num2str(),num2str(), num2str(obj.prior.alpha),num2str(obj.prior.kappa), obj.prior.m, num2str(obj.prior.v), num2str(obj.prior.Ws)};
                options.Interpreter = 'tex';
                if nargin > 1 && ishandle(src)
                    [result,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);
                    if ~Cancelled
                        obj.prior.alpha = str2double(result.alpha);
                        obj.prior.alphaA = str2double(result.alpha);
                        obj.prior.kappa = str2double(result.kappa);
                        obj.prior.m = result.m;
                        obj.prior.v = str2double(result.v);
                        obj.prior.W = str2double(result.Ws);
                        obj.prior.Ws = obj.prior.W;
                        obj.NTimes = str2double(result.NTimes);
                        obj.NOrder=str2double(result.NOrder);

                        
                        obj.maximumlikelihood = (result.List == 1);
                    end
                else
                    Cancelled = false;
                end                
                                        
                 if ~Cancelled 
                    priorSub.alpha = obj.prior.alpha;
                    priorSub.alphaA = obj.prior.alphaA;
                    priorSub.kappa = obj.prior.kappa;
                    priorSub.m = eval(obj.prior.m);
                    priorSub.v = obj.prior.v;
                    priorSub.M =  obj.prior.W;
                     
                    [labelposterior, modelHMM, obj.evidence, modelGMM,labellikelihood] = bayesHMM(dataTrans, obj.NOrder,obj.NTimes,priorSub);
                    if obj.maximumlikelihood
                        label = labellikelihood;
                    else
                        label = labelposterior;
                    end
                   
                    obj.selectedModel=modelHMM;
                    a = find(obj.spotsIncluded);
                    dataIdx =~isnan(dataTrans(1,:,:)); 
                    tempM = zeros([size(dataTrans,2) size(dataTrans,3)]);
                    tempM(dataIdx) =  label;
                    obj.fittedTraces = zeros([size(obj.spotsIncluded,1) size(idxData,2)]);
                    obj.fittedTraces(a(indexInital),:) = tempM';
                    
                    h1=figure;
                    dataIdx = ~isnan(dataTrans);
                    X = reshape(dataTrans(dataIdx),[size(dataTrans,1) sum(sum(sum(dataIdx)))./size(dataTrans,1)]);
                    
                    h=subplot(1,2,2,'parent',h1)
                    gplotmatrix2(X',[],label,['c' 'b' 'm' 'g' 'r'],[],[],false,[],{'Background' 'Intensity'},[],h);
                    title([sprintf('Model means are') mat2str(round(modelHMM.m*100)/100)...
                    sprintf('\nThe occupancy is ') mat2str(round(sum(obj.selectedModel.R,1)./size(obj.selectedModel.R,1)*100)/100) ])
                    h=subplot(1,2,1,'parent',h1)
                    h2 = plot(mean(obj.evidence,2),'ob','parent',h)
                    kopt = size(modelHMM.R,2);
                    line([kopt kopt],[min(min(obj.evidence)) max(max(obj.evidence))],...
                      'Color','r','parent',h);
                    axis tight           
                    title([sprintf('Selected Model N=%d',kopt)])
                 end                
%             end
            obj.selectMenu(2);
            obj.changeFilterPosition;             
       end
       
       function [traces2_allfr, ttb] = getTraceMatrix(obj)
            idx = find(obj.getSelectedMenu,1,'first');
            if isempty(idx)
                idx=1;
                obj.selectMenu(1);
            end
            switch(idx)
                case {1}
                    [traces2_in] = obj.getBinary(true(size(obj.spotsIncluded)));
                case{2}
                    traces2_in = obj.fittedTraces(true(size(obj.spotsIncluded)),:);
            end
            
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
            
            traces2_allfr = traces2_allfr(obj.spotsIncluded,:);

            se1 = strel('line', max(1,obj.handels.mdSlider.Value), 0);
            se2 = strel('line', max(1,obj.handels.mgSlider.Value), 0);

            A1 = imclose(traces2_allfr,se1);
            traces2_allfr = imopen(A1,se2);
            obj.spotNrs = find(obj.spotsIncluded);
           
            if isempty(obj.intergrationTime)
                obj.intergrationTime=1;
            end
    
            
            if ischar(obj.intergrationTime)
                ws = load(obj.intergrationTime,'vid')
                ttb = ws.vid.ttb;
            elseif isnumeric(obj.intergrationTime)
                ttb=1000*[0:1/obj.intergrationTime:1/obj.intergrationTime*(size(traces2_allfr,2))]; % This is an array containing relative time stamps for all frames
            end
       end
       
       function copyRastergramSettings(obj,obj2)
            obj.rastergramStartframe = obj2.rastergramStartframe; % write the answer before ; mark
            obj.sortType = obj2.sortType ; %none=0;duration of first=1;duration of last=2;binding of first=3;
            obj.sortState = obj2.sortState; % state to sort detult is 1
            obj.align = obj2.align; %none=0;first binding event=1;last departure event=2;
            obj.delete = obj2.delete; %true = 1 and false = 0
            obj.rastergramStartframe = obj2.rastergramStartframe; %frame to start rastergram
            obj.colorArray = obj2.colorArray; %rgb values of the color of each state e.g [1 1 1;0 0 1;1 0 0;1 0 0;0 1 0;1 0.4 0.4;0 0 0]
            obj.numberOfSpotsForRastergram = obj2.numberOfSpotsForRastergram; % empty for using all states
       end
       
        function setRastergramSettings(obj,src,~)
            Title = 'Rastergram Settings';

            %%%% SETTING DIALOG OPTIONS
            Options.Resize = 'on';
            Options.Interpreter = 'tex';
            Options.CancelButton = 'on';
            Options.ApplyButton = 'off';
            
            Option.Dim = 1; % Horizontal dimension in fields
            Prompt = {};
            Formats = {};

            Prompt(1,:) = {'Startframe Rastergrams','rastergramStartframe',[]};
            Formats(1,1).type = 'edit';
            Formats(1,1).format = 'text';
            Formats(1,1).size = [-1 0];
            Formats(1,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.rastergramStartframe = mat2str(max(1,obj.rastergramStartframe));
                       
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
            DefAns.delete = logical(obj.delete);

            Prompt(5,:) = {'Number of Spots','numberOfSpotsForRastergram',[]};
            Formats(5,1).type = 'edit';
            Formats(5,1).format = 'text';
            Formats(5,1).size = [-1 0];
            Formats(5,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.numberOfSpotsForRastergram = num2str(obj.numberOfSpotsForRastergram);  
            idx = find(obj.getSelectedMenu,1,'first');
            menuOneCheck = ~(obj.status == 1) || idx == 1;
            if ~menuOneCheck
                Prompt(6,:) = {'Sort & Align State','sortState',[]};
                Formats(6,1).type = 'edit';
                Formats(6,1).format = 'integer';
                Formats(6,1).size = [-1 0];
                Formats(6,1).span = [1 1];  % item is 1 field x 3 fields
                DefAns.sortState = obj.sortState;
            end
            
            Prompt(end+1,:) = {'Color Array','colorArray',[]};
            Formats(end+1,1).type = 'edit';
            Formats(end,1).format = 'text';
            Formats(end,1).size = [-1 0];
            Formats(end,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.colorArray = mat2str(max(0,min(1,obj.colorArray))); %[pwd '/Dark.tif'];


            [result,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);

            if ~Cancelled
                obj.sortType = result.List-1;
                obj.align = result.align-1;
                obj.delete = result.delete;
                obj.rastergramStartframe = eval(result.rastergramStartframe);
                obj.colorArray = eval(result.colorArray);
                if menuOneCheck
                    obj.sortState = 1;
                else
                    obj.sortState = result.sortState;
                end
                if isempty(result.numberOfSpotsForRastergram)
                    obj.numberOfSpotsForRastergram =[];
                else
                    obj.numberOfSpotsForRastergram = str2num(result.numberOfSpotsForRastergram);
                end
            end
        end

        function analyze(obj,src,~)
            idx = find(obj.getSelectedMenu,1,'first');
            if ~(obj.status == 1) || idx == 1
                analyze@guiAnalyzeDynamics(obj);  
            end
            obj.selectMenu(1);
            obj.sortState=1;
            obj.changeFilterPosition;              
        end

        function plotImages(obj,~,~)
            plotImages@guiAnalyzeDynamics(obj);            
            obj.changeFilterPosition;
        end
        
        function loadvars(obj,src,updateplot)
            if nargin < 3 || isempty(updateplot) || isobject(updateplot)
                updateplot=true;
            end
            updateplot=updateplot&(~isempty(obj.dataObject) && (isfield(obj.dataObject,'status') && obj.dataObject.status));
            obj.draw3Dmovie=updateplot;
            
            ws = loadvars@guiAnalyzeDynamics(obj,src,updateplot);
            gabDurvec(1) = ws.mdSliderValue;
            gabDurvec(2) = ws.mgSliderValue;
            obj.setGabDurSliders(gabDurvec);
        end
        
        function savevars(obj,src,~)
            vs = savevars@guiAnalyzeDynamics(obj);
            gabDurvec = obj.getGabDurSliders;
            vs.mdSliderValue = gabDurvec(1);
            vs.mgSliderValue = gabDurvec(2);
            
            if nargin > 1 && isobject(src)
                [file, pathname] = uiputfile([vs.class '.mat'] );
                if file ~= 0
                    answer{1} = fullfile(pathname, file);
                else
                    answer=[];
                end
                if ~isempty(answer)
                    save(answer{1},'vs','-v7.3')
                end
            elseif ischar(src)
                save(src,'vs','-v7.3')
            end
        end
        
        function setGabDurSliders(obj,gabDurVec)
             obj.handels.mdSlider.Value = round(max(1,gabDurVec(1)));
             obj.handels.mgSlider.Value = round(max(1,gabDurVec(2)));
             obj.updateFilter;
        end
        
        function gabDurVec = getGabDurSliders(obj)
            gabDurVec(1) = round(max(1,obj.handels.mdSlider.Value));
            gabDurVec(2) = round(max(1,obj.handels.mgSlider.Value));
        end
        
        
        function exportData(obj,src,~)
            traces2_allfr = obj.getTraceMatrix;
            
            [file, pathname] = uiputfile([class '.mat'] );
            if file ~= 0
                answer{1} = fullfile(pathname, file);
            else
                answer=[];
            end
            if ~isempty(answer)
                data.Trace = traces2_allfr;
                data.Background = obj.bg1I(obj.spotsIncluded,:,2-0)*(obj.PSFSigma*2)^2;
                data.Intensity = obj.photons1I(obj.spotsIncluded,:,2-0);
                data.DxSq =  obj.deltay(obj.spotsIncluded,:).^2+obj.deltax(obj.spotsIncluded,:).^2;
                data.Sigma =  obj.Sigma1I(obj.spotsIncluded,:,2-0);
                data.spotIDs = obj.spotNrs;
                save(answer{1},'data','-v7.3')
            end
        end
    end
end