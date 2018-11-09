classdef guiDataCollector < handle
   properties
        objects
        class = 'guiDataCollector' 
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
       
        function obj = guiDataCollector(varargin)
            if isa(varargin{1},'guiDataLoader')
                uimenu(varargin{1}.handels.parentMenu,'Label','Load All','Separator','on','Callback',@obj.loadvars);
                uimenu(varargin{1}.handels.parentMenu,'Label','Save  All','Callback',@obj.savevars);
            else
                error('First argument should be a guiDataLoader')
            end
            
            obj.objects = varargin;
        end
        
        function savevars(obj,src,~)
            class = obj.class;
            if ishandle(src)
                [file, pathname] = uiputfile([obj.class '.mat']);
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
                dlg = msgbox(['Saving Data files']);
                for i=2:length(obj.objects)
                    vs{i} = obj.objects{i}.savevars;
                end           

                save(answer{1},'class','vs','-v7.3')

                obj.objects{1}.savevars([answer{1} ]);

                if ishghandle(dlg)
                    delete(dlg);
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
                    obj.objects{1}.loadvars(['guiDataLoader' answer{1} ]); 
                    ws = load(answer{1});
                    if ~strcmpi(ws.class,obj.class)
                        error('This is no guiSpotDetector class')
                    else
                        for i=2:length(obj.objects)
                            obj.objects{i}.loadvars(ws.vs{i},updateplot);
                        end
                    end
            end
       end
   end
end