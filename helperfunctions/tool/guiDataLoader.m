classdef guiDataLoader < handle
   properties 
    Figure
    Axes
    CLimit
    
    dataCrop
    
    NFramesGroup=1
    StartFrameSummation=1
    
    name
    hsl1
    s
    htextbox
    filename
    C
    
    numberOfChannels
    dataCam
    
    beadDataCam1
    beadDataCam2
    
    img1
    
    sourceFileLoad
    
    %results
    outCam1
    outCam2
    bgMean
    
    getCount
    
    handels
    status
    class = 'guiDataLoader'
    loadedCal = 0
    loadedExp = 0
    reloadData = 1
    reloadCal = 1
    reloadalgn =1
    private

    fileStructure
    flipCams

   end
   events
       DataChange
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
       
       function intT = getFramesPerSecond(obj,src,~)
           intT = obj.NFramesGroup/obj.fileStructure.intT;
       end

	function obj = selectExpPath(obj,src,~)
        obj.reloadData=1;
            if isstruct(src) && ~ishandle(src) 
                obj.fileStructure.Data = src.Data;

                if isfield(src,'lambdar')
                    obj.fileStructure.lambdar = src.lambdar;
                end
                if isfield(src,'lambdal')
                    obj.fileStructure.lambdal = src.lambdal;
                end
                if isfield(src,'flipCams')
                    obj.fileStructure.flipCams = src.flipCams;
                    obj.flipCams = src.flipCams;
                end
                if isfield(src,'intT')
                    obj.fileStructure.intT = src.intT;
                end
                obj.loadedExp=1;

                if obj.checkStatus == 1
                    obj.plotImages;
                end            
            else
            Title = 'Load files';

            %%%% SETTING DIALOG OPTIONS
            Options.Resize = 'on';
            Options.Interpreter = 'tex';
            Options.CancelButton = 'on';
            Options.ApplyButton = 'off';
            Option.Dim = 1; % Horizontal dimension in fields

            Prompt = {};
            Formats = {};
            % DefAns = struct([]);

            Prompt(1,:) = {'Data','Data',[]};
            Formats(1,1).type = 'edit';
            Formats(1,1).format = 'file';
            Formats(1,1).size = [-1 0];
            Formats(1,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.Data = pwd; 
            
            Prompt(2,:) = {'Left Wavelength [nm]','lambdal',[]};
            Formats(2,1).type = 'edit';
            Formats(2,1).format = 'integer';
            Formats(2,1).size = [-1 0];
            Formats(2,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.lambdal = 647; 

            Prompt(3,:) = {'Right Wavelength [nm]','lambdar',[]};
            Formats(3,1).type = 'edit';
            Formats(3,1).format = 'integer';
            Formats(3,1).size = [-1 0];
            Formats(3,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.lambdar = 550; 

            Prompt(4,:) = {'frames per second','intT',[]};
            Formats(4,1).type = 'edit';
            Formats(4,1).format = 'integer';
            Formats(4,1).size = [-1 0];
            Formats(4,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.intT = 1; 
            
            Prompt(5,:) = {'flipCams','flipCams',[]};
            Formats(5,1).type = 'check';
            Formats(5,1).format = 'logical';
            Formats(5,1).size = [-1 0];
            Formats(5,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.flipCams = false;
            
            Prompt(6,:) = {'Number of frames to group','NFramesGroup',[]};
            Formats(6,1).type = 'edit';
            Formats(6,1).format = 'integer';
            Formats(6,1).size = [-1 0];
            Formats(6,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.NFramesGroup = 1; 

            
            [ExpfileStructure,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);
            obj.fileStructure.Data = ExpfileStructure.Data;
            
            obj.fileStructure.flipCams = ExpfileStructure.flipCams;
            obj.fileStructure.intT = ExpfileStructure.intT;
            obj.fileStructure.lambdal = ExpfileStructure.lambdal;
            obj.fileStructure.lambdar = ExpfileStructure.lambdar;
            
            obj.flipCams = obj.fileStructure.flipCams;
            if ~Cancelled
                obj.loadedExp=1;
            end
            
            if obj.checkStatus
                obj.loadImage;
            end
            end
            
    end
       
	function obj = selectCalPath(obj,src,~)
                obj.reloadCal=1;
                obj.reloadalgn = 1;
        if isstruct(src) && ~ishandle(src)  
            obj.fileStructure.Grid = src.Grid;
            obj.fileStructure.Dark =  src.Dark;
            obj.fileStructure.BeadData =  src.BeadData;
            if isfield(src,'lambdar')
                obj.fileStructure.lambdar = src.lambdar;
            end
            if isfield(src,'lambdal')
                obj.fileStructure.lambdal = src.lambdal;
            end
            if isfield(src,'flipCams')
                obj.fileStructure.flipCams = src.flipCams;
                obj.flipCams = src.flipCams;
            end
            if isfield(src,'intT')
                obj.fileStructure.intT = src.intT;
            end
            if ~isfield(src,'rg') || isempty(src.rg)
                obj.fileStructure.rg = -1;
            elseif isa(src.rg,'char') && strcmpi(src.rg,'[0 Inf]')
                obj.fileStructure.rg = eval(obj.fileStructure.rg);    
            elseif isa(src.rg,'double')
                obj.fileStructure.rg = obj.fileStructure.rg;
            end

            obj.loadedCal=1;

            if obj.checkStatus
                obj.loadImage;
            end   
        else
            Title = 'Load files';

            %%%% SETTING DIALOG OPTIONS
            Options.Resize = 'on';
            Options.Interpreter = 'tex';
            Options.CancelButton = 'on';
            Options.ApplyButton = 'off';
            Option.Dim = 1; % Horizontal dimension in fields

            Prompt = {};
            Formats = {};

            Prompt(1,:) = {'Grid','Grid',[]};
            Formats(1,1).type = 'edit';
            Formats(1,1).format = 'file';
            Formats(1,1).size = [-1 0];
            Formats(1,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.Grid = pwd;

            Prompt(2,:) = {'Dark','Dark',[]};
            Formats(2,1).type = 'edit';
            Formats(2,1).format = 'file';
            Formats(2,1).size = [-1 0];
            Formats(2,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.Dark = pwd;
            
            Prompt(3,:) = {'Intensity Range','rg',[]};
            Formats(3,1).type = 'edit';
            Formats(3,1).format = 'text';
            Formats(3,1).size = [-1 0];
            Formats(3,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.rg = '[0 Inf]'; 


            Prompt(4,:) = {'Bead','BeadData',[]};
            Formats(4,1).type = 'edit';
            Formats(4,1).format = 'file';
            Formats(4,1).size = [-1 0];
            Formats(4,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.BeadData =  pwd;

         

            Prompt(5,:) = {'flipCams','flipCams',[]};
            Formats(5,1).type = 'check';
            Formats(5,1).format = 'logical';
            Formats(5,1).size = [-1 0];
            Formats(5,1).span = [1 1];  % item is 1 field x 3 fields
            DefAns.flipCams = false;

            [calfileStructure,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);

            obj.fileStructure.Grid = calfileStructure.Grid;
            obj.fileStructure.Dark = calfileStructure.Dark;
            obj.fileStructure.BeadData = calfileStructure.BeadData;
            obj.fileStructure.flipCams = calfileStructure.flipCams;
            if strcmpi(calfileStructure.rg,'[0 Inf]')
                obj.fileStructure.rg = -1;
            else
                obj.fileStructure.rg = eval(calfileStructure.rg);
            end

            if ~Cancelled
                obj.loadedCal =1;
            end
            if obj.checkStatus == 1
                obj.plotImages;
            end
        end
       end
 
       
       function obj = selectPath(obj,src,~)
            if isstruct(src) && ~ishandle(src)
               obj.clear;
               obj.fileStructure = src;    
               obj.flipCams = obj.fileStructure.flipCams;
               obj.loadedExp=1;
               obj.loadedCal=1;
               if obj.checkStatus == 1
                    obj.loadImage;
               end   
            else
                Title = 'Load files';

                %%%% SETTING DIALOG OPTIONS
                Options.Resize = 'on';
                Options.Interpreter = 'tex';
                Options.CancelButton = 'on';
                Options.ApplyButton = 'off';
                Option.Dim = 1; % Horizontal dimension in fields

                Prompt = {};
                Formats = {};
                
                Prompt(1,:) = {'Grid','Grid',[]};
                Formats(1,1).type = 'edit';
                Formats(1,1).format = 'file';
                Formats(1,1).size = [-1 0];
                Formats(1,1).span = [1 1];  % item is 1 field x 3 fields
                DefAns.Grid = pwd;

                Prompt(2,:) = {'Dark','Dark',[]};
                Formats(2,1).type = 'edit';
                Formats(2,1).format = 'file';
                Formats(2,1).size = [-1 0];
                Formats(2,1).span = [1 1];  % item is 1 field x 3 fields
                DefAns.Dark = pwd;                
                    
                Prompt(3,:) = {'Intensity Range','rg',[]};
                Formats(3,1).type = 'edit';
                Formats(3,1).format = 'text';
                Formats(3,1).size = [-1 0];
                Formats(3,1).span = [1 1];  % item is 1 field x 3 fields
                DefAns.rg = '[0 Inf]'; 



                Prompt(4,:) = {'Bead','BeadData',[]};
                Formats(4,1).type = 'edit';
                Formats(4,1).format = 'file';
                Formats(4,1).size = [-1 0];
                Formats(4,1).span = [1 1];  % item is 1 field x 3 fields
                DefAns.BeadData =  pwd;

                Prompt(5,:) = {'Data','Data',[]};
                Formats(5,1).type = 'edit';
                Formats(5,1).format = 'file';
                Formats(5,1).size = [-1 0];
                Formats(5,1).span = [1 1];
                DefAns.Data = pwd;
                
                Prompt(6,:) = {'Number of frames to group','NFramesGroup',[]};
                Formats(6,1).type = 'edit';
                Formats(6,1).format = 'integer';
                Formats(6,1).size = [-1 0];
                Formats(6,1).span = [1 1]; 
                DefAns.NFramesGroup = 1; 

                Prompt(7,:) = {'Left Wavelength [nm]','lambdal',[]};
                Formats(7,1).type = 'edit';
                Formats(7,1).format = 'integer';
                Formats(7,1).size = [-1 0];
                Formats(7,1).span = [1 1];
                DefAns.lambdal = 647;

                Prompt(8,:) = {'Right Wavelength [nm]','lambdar',[]};
                Formats(8,1).type = 'edit';
                Formats(8,1).format = 'integer';
                Formats(8,1).size = [-1 0];
                Formats(8,1).span = [1 1];
                DefAns.lambdar = 550;

                Prompt(9,:) = {'frames per second','intT',[]};
                Formats(9,1).type = 'edit';
                Formats(9,1).format = 'integer';
                Formats(9,1).size = [-1 0];
                Formats(9,1).span = [1 1];
                DefAns.intT = 1; 

                Prompt(10,:) = {'flipCams','flipCams',[]};
                Formats(10,1).type = 'check';
                Formats(10,1).format = 'logical';
                Formats(10,1).size = [-1 0];
                Formats(10,1).span = [1 1];
                DefAns.flipCams = false;

                [fileStructureInp,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);
         
                if ~Cancelled
                    obj.clear;
                    obj.fileStructure=fileStructureInp;
                    obj.NFramesGroup=fileStructureInp.NFramesGroup;
                    obj.flipCams = obj.fileStructure.flipCams;
                    obj.loadedCal=1;
                    obj.loadedExp=1;
                    
                    if strcmpi(fileStructureInp.rg,'[0 Inf]')
                        obj.fileStructure.rg = -1;
                    else
                        obj.fileStructure.rg = eval(fileStructureInp.rg);
                    end
                end
                
                if obj.checkStatus == 1
                    obj.loadImage;
                end
            end
       end
       
      function obj = guiDataLoader(h,fileStructure,flipCams,NFramesGroup,StartFrameSummation)
          obj.dataCam{1}=[];
          obj.dataCam{2}=[];
          obj.dataCam{3}=[];
          obj.dataCam{4}=[];
          obj.status=0;
          
          if nargin > 2
              obj.NFramesGroup = NFramesGroup;
              obj.StartFrameSummation=StartFrameSummation;
          end
          
          if nargin > 1 && isfield(fileStructure,'NFramesGroup')
              obj.NFramesGroup = fileStructure.NFramesGroup;
          end
          
          if nargin > 1 && isfield(fileStructure,'StartFrameSummation')
              obj.StartFrameSummation = fileStructure.StartFrameSummation;
          end
          
          obj.getCount=zeros(4,1);
            obj.numberOfChannels = 2;
            
            htab1 = uitab(h, 'Title', 'Data');
            obj.Figure = htab1;
            obj.dataCrop = [];
            
            obj.handels.parentMenu = uimenu('Label','Data');
            uimenu(obj.handels.parentMenu,'Label','Clear','Callback',@obj.clear);
            uimenu(obj.handels.parentMenu,'Label','Load Cal Data','Separator','on','Callback',@obj.selectCalPath);
            uimenu(obj.handels.parentMenu,'Label','Load Exp Data','Callback',@obj.selectExpPath);          
            
            if nargin > 1 && ~isempty(fileStructure)
               obj.fileStructure = fileStructure;
            else
                uimenu(obj.handels.parentMenu,'Label','Load All Data','Callback',@obj.selectPath);
            end
            
            if (nargin > 2 && ~isempty(flipCams)) 
                obj.flipCams = flipCams;
            elseif (nargin > 1 && ~isempty(fileStructure) && isfield(fileStructure,'flipCams'))
                obj.flipCams=fileStructure.flipCams;
            else
                obj.flipCams=false;
            end
            
            uimenu(obj.handels.parentMenu,'Label','Load Callibration','Separator','on','Callback',@obj.loadcal);
            uimenu(obj.handels.parentMenu,'Label','Save Callibration','Callback',@obj.savecal);

            uimenu(obj.handels.parentMenu,'Label','Load Data & Callibration','Separator','on','Callback',@obj.loadvars);
            uimenu(obj.handels.parentMenu,'Label','Save Data & Callibration','Callback',@obj.savevars);

            slmin = 1;
            slmax = 100;
            obj.hsl1 = uicontrol('Parent', obj.Figure,'Style','slider','Min',slmin,'Max',slmax,...
                            'SliderStep',[1 1]./(slmax-slmin),'Value',51,...
                            'Units','normalized','Position',[ 0.5 0.01 0.4 0.05]);
                        
            slmin = 0;
            slmax = 100;
            obj.s = uicontrol('Parent', obj.Figure,'Style','Slider',...
            'Units','normalized','Position',[0.95 0 0.05 1],'Min',slmin,'Max',slmax,...
                  'Value',51);                         
            
            set(obj.s,'Callback',@obj.plotImages)
            set(obj.hsl1,'Callback',@obj.plotImages)
            
                      obj.htextbox = uicontrol(... 
                    'Parent', obj.Figure,...
                    'Units','normalized',...
                    'Position',[.01,-.01,.45,1],...
                    'BackgroundColor','w',...
                    'HorizontalAlignment','left',...
                    'FontSize',12,...
                    'Style','text',...
                    'String',[ '']);
            
            if nargin > 1 && ~isempty(fileStructure)
                 obj.loadImage;
            end
      end
      function loadImage(obj,src,~)

        if obj.reloadCal || isempty(obj.outCam1) || isempty(obj.outCam2) || isempty(obj.bgMean)
            obj.loadImageCallibration;
            obj.reloadData=0;
        end
        if  obj.reloadalgn || isempty(obj.dataCam{1}) ||  isempty(obj.dataCam{2}) 
            obj.loadAlignmentData;
        end
        if  obj.reloadData || isempty(obj.dataCam{3}) ||  isempty(obj.dataCam{3}) 
            obj.loadData;                       	
        end
        if obj.checkStatus == 1
            obj.plotImages;
        end
      end
      
      function st =  checkStatus(obj,src,~)
          if obj.loadedCal && obj.loadedExp
              obj.status=1;
          else
              obj.status=0.5;
          end
          st = obj.status;
      end
         
        function loadImageCallibration(obj,src,~)

            %% Callibration Data
            a = bfopen( [ obj.fileStructure.Grid]);
            geller = double(cell2mat(permute(a{1}(1:min(100,end),1),[3 2 1])));
            
            a = bfopen( [ obj.fileStructure.Dark] );
            bg = double(cell2mat(permute(a{1}(1:min(100,end),1),[3 2 1])));
            
            gellerCam2 = geller(1:512,1:512,:);
            gellerCam1 = geller(1:512,513:1024,:);
            clear geller
            
            bgCam2 = bg(1:512,1:512,:);
            bgCam1 = bg(1:512,513:1024,:);
            obj.bgMean = mean(bg,3);
            clear bg
            
            obj.outCam1 = cal_readnoise(gellerCam1, bgCam1,100,obj.fileStructure.rg);
            obj.outCam2 = cal_readnoise(gellerCam2, bgCam2,100,obj.fileStructure.rg);
            obj.loadedCal = 1;
        end
        function color = getColor(obj,idx)
            if mod(idx,2)==0 %  iseven
                color = spectrumRGB(obj.fileStructure.lambdal);
            elseif mod(idx,2)==1 % isodd
                color = spectrumRGB(obj.fileStructure.lambdar);
            else
                error('Channel not odd or even!')
            end
        end
        function data = getData(obj,idx,pixelx,pixely,pixelz)
            
            obj.getCount(idx)=obj.getCount(idx)+1;
            if obj.checkStatus
                if isempty(obj.dataCam{idx})
                    variable{1} = 'dataCam1';
                    variable{2} = 'dataCam2';
                    variable{3} = 'beadDataCam1';
                    variable{4} = 'beadDataCam2';
                    dlg = msgbox(['Loading ' variable{idx} ' operation in progress...']);
                    
                    if ~isempty(obj.sourceFileLoad)
                        ws = load(obj.sourceFileLoad,variable{idx});
                        eval(['obj.dataCam{idx} = ws.' variable{idx} ';']);  
                    else
                        if idx == 1 || idx == 2
                            obj.loadData;                       	
                        elseif idx == 3 || idx == 4
                            obj.loadAlignmentData;
                        else      
                            msgbox(['Error idx - no file to load data from!']);
                        end
                    end
                   
                    if ishghandle(dlg)
                        delete(dlg);
                    end 
                end

                if nargin < 3 || isempty(pixelx)
                     pixelx = 1:size(obj.dataCam{idx},1);
                end
                if nargin < 3 || isempty(pixely)
                    pixely = 1:size(obj.dataCam{idx},2);
                end
                if nargin < 3 || isempty(pixelz)
                    pixelz = 1:size(obj.dataCam{idx},3);
                end

                data = obj.dataCam{idx}(pixelx,pixely,pixelz);
            else
                data = [];
            end
        end
        
        function loadAlignmentData(obj,src,~)
            nFrames = 10;

            a = bfopen( [ obj.fileStructure.BeadData] );
            data = double(cell2mat(permute(a{1}(1:min(end,nFrames),1),[3 2 1])));
            
            data=(data-repmat(obj.bgMean,[1 1 size(data,3)]));
         
            obj.dataCam{4-obj.flipCams} = data(1:512,1:512,1:min(end,nFrames))*obj.outCam1(2); % swapped Cam1 and Cam2 - VS
            obj.dataCam{3+obj.flipCams} = data(1:512,513:1024,1:min(end,nFrames))*obj.outCam2(2);   
            obj.loadedCal =1;
        end
         
        function loadData(obj,src,~)
 
            nFrames = 1e6; 
            if isfield(obj.fileStructure,'Data')
                a = bfopen( [ obj.fileStructure.Data] );
                data = double(cell2mat(permute(a{1}(1:min(end,nFrames),1),[3 2 1])));
                clear a
                data=(data-repmat(obj.bgMean,[1 1 size(data,3)]));

                obj.dataCam{2-obj.flipCams} = data(1:512,1:512,1:min(end,nFrames))*obj.outCam1(2); % swapped Cam1 and Cam2 - VS
                obj.dataCam{1+obj.flipCams} = data(1:512,513:1024,1:min(end,nFrames))*obj.outCam2(2);
                clear data

                if obj.flipCams
                    temp = obj.fileStructure.lambdal;
                    obj.fileStructure.lambdal = obj.fileStructure.lambdar;
                    obj.fileStructure.lambdar = temp;
                    clear temp
                end

                if obj.NFramesGroup > 1 
                        if obj.StartFrameSummation > 1
                            endPoint=floor(size(obj.dataCam{2}(:,:,obj.StartFrameSummation:end),3)/obj.NFramesGroup)*obj.NFramesGroup;
                            idxSummed = obj.StartFrameSummation:obj.StartFrameSummation+endPoint-1;
                            B = reshape(permute(obj.dataCam{1}(:,:,idxSummed),[3 1 2 ]),[obj.NFramesGroup length(idxSummed)/obj.NFramesGroup size(obj.dataCam{1},1) size(obj.dataCam{1},2)]);
                            obj.dataCam{1} = cat(3,obj.dataCam{1}(:,:,1:obj.StartFrameSummation-1),permute(squeeze(sum(B,1)),[2 3 1]));

                            B = reshape(permute(obj.dataCam{2}(:,:,idxSummed),[3 1 2 ]),[obj.NFramesGroup length(idxSummed)/obj.NFramesGroup size(obj.dataCam{1},1) size(obj.dataCam{1},2)]);
                            obj.dataCam{2} = cat(3,obj.dataCam{2}(:,:,1:obj.StartFrameSummation-1),permute(squeeze(sum(B,1)),[2 3 1]));
                        else
                            B = reshape(permute(obj.dataCam{1}(:,:,1:floor(size(obj.dataCam{2},3)/obj.NFramesGroup)*obj.NFramesGroup),[3 1 2 ]),[obj.NFramesGroup floor(size(obj.dataCam{2},3)/obj.NFramesGroup) size(obj.dataCam{1},1) size(obj.dataCam{1},2)]);
                            obj.dataCam{1} = permute(squeeze(sum(B,1)),[2 3 1]);

                            B = reshape(permute(obj.dataCam{2}(:,:,1:floor(size(obj.dataCam{2},3)/obj.NFramesGroup)*obj.NFramesGroup),[3 1 2 ]),[obj.NFramesGroup floor(size(obj.dataCam{2},3)/obj.NFramesGroup) size(obj.dataCam{1},1) size(obj.dataCam{1},2)]);
                            obj.dataCam{2} = permute(squeeze(sum(B,1)),[2 3 1]);
                        end
                     if rem(size(obj.dataCam{2},3)/obj.NFramesGroup,1) ~= 0
                        warning('total of frams should be divisible by frame group number')
                     end
                end
                maxNumberOfImages=4-1;
                set(obj.hsl1, 'SliderStep', [1/maxNumberOfImages , 10/maxNumberOfImages ]);       
                obj.loadedExp=1;
                notify(obj,'DataChange');
            else
                disp('No data to load; set fileStructure')
            end
        end


          function plotImages(obj,src,~)
            set(0,'CurrentFigure',getParentFigure(obj.Figure));
            channelPos = round(get(obj.hsl1,'Value'));
            chan = round((obj.numberOfChannels*2-1)/100*channelPos)+1;
            if ~isempty(squeeze(obj.getData(chan)))
                    switch chan
                        case 1
                            str = 'dataCam{1}';
                            name = 'Exp Camera 1';
                        case 2
                            str = 'dataCam{2}';
                            name = 'Exp Camera 2';
                        case 3
                            str = 'dataCam{3}';
                            name = 'Reg Data Camera 1';
                        case 4
                            str = 'dataCam{4}';
                            name = 'Reg Data Camera 2';
                        otherwise 
                            disp('?')
                    end
                        
                    
                    sliderPos = round(get(obj.s,'Value'));
                    slice = round((size(eval(['obj.' str]),3)-1)/100*sliderPos);

                    x1=0;
                    y1=0;

                    hAxis2 = axes('Parent', obj.Figure,'Position', [x1+0.501 y1 0.45 1]);
                    axis off; % Turn off tick marks, etc.
                    gimr = imadjust(mat2gray(double(stretch(eval(['obj.getData(' num2str(chan) ',[],[],slice+1)']))))); %flip(flip(A(:,:,ii+1),1),2));
                    
                    color = obj.getColor(chan);
                    obj.img1 = imshow(cat(3,color(1).*gimr,color(2).*gimr,color(3).*gimr));
                     
                    rs = obj.getSummary;
                    txt1 = ['On display = ' name ', slice = ' num2str(slice)]; 
                    
                    set(obj.htextbox,'String',{txt1;rs{1}; rs{2}; rs{3}; rs{4}});
            end    
          end
          
          function rs = getSummary(obj)
                [pathstr,name,ext] = fileparts(obj.fileStructure.Grid);
                txt1 = ['Grid Filename = ' name];
 
                [pathstr,name,ext] = fileparts(obj.fileStructure.Dark);
                txt2 = ['Dark Filename = ' name];

                [pathstr,name,ext] = fileparts(obj.fileStructure.BeadData);
                txt3 = ['Bead Filename = ' name];

                [pathstr,name,ext] = fileparts(obj.fileStructure.Data);
                txt4 = ['Data Filename = ' name];

                rs = {txt1;txt2;txt3;txt4};
          end
          
          function clear(obj,src,~)
            obj.dataCam{1}=[];
            obj.dataCam{2}=[];

            obj.dataCam{3}=[];
            obj.dataCam{4}=[];
            
            obj.status=0;
            obj.loadedCal =0;
            obj.loadedExp=0;
            obj.flipCams =0;
            obj.fileStructure=[];
            set(obj.htextbox,'String','');
            
            delete(obj.img1)
            notify(obj,'DataChange');
          end
          
          function savevars(obj,src,~)
               
            class = obj.class;       
            dataCam1 = obj.dataCam{1};           
            dataCam2 = obj.dataCam{2};           
            beadDataCam1 = obj.dataCam{3};           
            beadDataCam2 = obj.dataCam{4}; 
            fileStructure=obj.fileStructure;
            
            if nargin > 1 && isobject(src)
                    [file, pathname] = uigetfile([obj.class '.mat'] );
                    if file ~= 0
                        answer{1} = fullfile(pathname, file);
                    else
                        answer=[];
                    end
                    if ~isempty(answer)
                        save(answer{1},'class','dataCam1','dataCam2','beadDataCam1','beadDataCam2','fileStructure','-v7.3')
                    end
            elseif ischar(src)
                   save(src,'class','dataCam1','dataCam2','beadDataCam1','beadDataCam2','fileStructure','-v7.3')
            end
        end

        function loadvars(obj,src,~)
            if ishandle(src)
                [file, pathname] = uigetfile([obj.class '.mat'] );
                if file ~= 0
                    answer{1} = fullfile(pathname, file);
                else
                    answer=[];
                end
                    
            elseif ischar(src)
                answer{1} = src;
            end
            if ~isempty(answer)
                ws = load(answer{1},'class','fileStructure');
                obj.sourceFileLoad=answer{1};

                class = obj.class;
                if ~strcmpi(ws.class,class)
                    error(['This is no ' class ' class'])
                else
                    obj.fileStructure = ws.fileStructure;

                    rs = obj.getSummary;

                    set(obj.htextbox,'String',{rs{1}; rs{2}; rs{3}; rs{4}});
                end
                obj.loadedCal =1;
                obj.loadedExp=1;
            end
        end
        
        function savecal(obj,src,~)
               
            class = obj.class;                 
            beadDataCam1 = obj.dataCam{3};           
            beadDataCam2 = obj.dataCam{4};
            CalfileStructure  = obj.fileStructure;
            CalfileStructure.Data = [];
            
            bgMean  = obj.bgMean;            
            outCam1 = obj.outCam1;
            outCam2 = obj.outCam2;

            if nargin > 1 && isobject(src)
                    [file, pathname] = uiputfile([obj.class '.mat'] );
                    if file ~= 0
                        answer{1} = fullfile(pathname, file);
                    else
                        answer=[];
                    end
                    if ~isempty(answer)
                        save(answer{1},'class','beadDataCam1','beadDataCam2','CalfileStructure','bgMean','outCam1','outCam2','-v7.3')
                    end
            elseif ischar(src)
                   save(src,'class','-v7.3')
            end
        end

        function loadcal(obj,src,~)
            if ishandle(src)
                [file, pathname] = uigetfile([obj.class '.mat'] );
                if file ~= 0
                    answer{1} = fullfile(pathname, file);
                else
                    answer=[];
                end
            elseif ischar(src)
                answer{1} = src;
            end
            if ~isempty(answer)            
                ws = load(answer{1});
    
                class = obj.class;
                if ~strcmpi(ws.class,class)
                    error(['This is no ' class ' class'])
                else
           
                    obj.dataCam{3} = ws.beadDataCam1;     
                    obj.dataCam{4} = ws.beadDataCam2;

                    obj.fileStructure.Grid = ws.CalfileStructure.Grid;
                    obj.fileStructure.Dark =  ws.CalfileStructure.Dark;
                    obj.fileStructure.BeadData =  ws.CalfileStructure.BeadData;
                    obj.fileStructure.flipCams = ws.CalfileStructure.flipCams;
                    obj.fileStructure.lambdar = ws.CalfileStructure.lambdar;
                    obj.fileStructure.lambdal = ws.CalfileStructure.lambdal;

                    obj.bgMean = ws.bgMean;            
                    obj.outCam1 = ws.outCam1;
                    obj.outCam2 =  ws.outCam2;
                    obj.loadedCal =1;
            
                end
            end
        end
   end
end