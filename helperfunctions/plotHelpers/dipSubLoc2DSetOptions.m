function options = dipSubLoc2DSetOptions(varargin)
% dipSubLoc2DSetOptions   set options for dipSubLoc2D
% 
% options = dipSubLoc2DSetOptions;
% Creates a structure of options for dipSubLoc2D set to default values.
%
% options = dipSubLoc2DSetOptions('param1',value1,'param2',value1,...)
% Creates a structure of options for dipSubLoc2D seting named parameters to
% specified values and all other parameters to their default.
% 
% options = dipSubLoc2DSetOptions(oldopts,'param1',value1,'param2',value2,...)
% Creates a copy of oldopts with the named parameters altered.

options = struct;

if (nargin == 0 || ischar(varargin{1}))
    options = struct('im',[],...
        'BoxCenters',[],...
        'BoxSize',[],...
        'BoxColor',[],...
        'plotBoxes',1,...
        'fitCoords',[],...
        'fitCoordsMarker',1,...
        'fitColor',[],...
        'h',[],...
        'linewidth',2,...
        'markersize',[5 4],...
        'updateImage',1,...
        'marker','o');
end
if(nargin > 0)
    if(ischar(varargin{1}))
        i = 1;
        while (i < length(varargin))
            options = setfield(options,varargin{i},varargin{i+1});
            i = i+2;
        end
    else if(isstruct(varargin{1}))
            if(isstruct(varargin{1}))
                options = varargin{1};
                j = 2;
                while (j < length(varargin))
                    options = setfield(options,varargin{j},varargin{j+1});
                    j = j+2;
                end
            else disp('error - in calling function');
            end
        end
    end
end


return