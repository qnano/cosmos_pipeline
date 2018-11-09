function H=shadedErrorBar(x,y,errBar,lineProps,transparent)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Error checking    
error(nargchk(3,5,nargin))


%Process y using function handles if needed to make the error bar
%dynamically
if iscell(errBar) && ~isvector(y)
    fun1=errBar{1};
    fun2=errBar{2};
    errBar=fun2(y);
    y=fun1(y);
elseif ~iscell(errBar) && isvector(y)
    y=y(:)';
else
    error('2nd and 3rd input arguments are not compatible')
end

if isempty(x)
    x=1:length(y);
else
    x=x(:)';
end

if length(x) ~= length(y)
    error('inputs x and y are not of equal lengths')
end


%If only one error bar is specified then we will mirror it, turning it into
%both upper and lower bars. 
if length(errBar)==length(errBar(:))
    errBar=repmat(errBar(:)',2,1);
else
    f=find(size(errBar)==2);
    if isempty(f), error('errBar has the wrong size'), end
    if f==2, errBar=errBar'; end
end

if length(x) ~= length(errBar)
    error('inputs x and y must have the same length as errBar')
end


%Set default options
defaultProps={'-k'};
if nargin<4 || isempty(lineProps)
    lineProps=defaultProps; 
end
if ~iscell(lineProps)
    lineProps={lineProps}; 
end


if nargin<5 || ~isnumeric(transparent)
    transparent=0; 
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Plot the main line. We plot this first in order to extract the RGB values
% for the line colour. I am not aware of a function that does this.
H.mainLine=plot(x,y,lineProps{:});


% Work out the color of the shaded region and associated lines
% Using alpha requires the render to be openGL and so you can't
% save a vector image. On the other hand, you need alpha if you're
% overlaying lines. We therefore provide the option of choosing alpha 
% or a de-saturated solid colour for the patch surface.

col=get(H.mainLine,'color');
edgeColor=col+(1-col)*0.55;
patchSaturation=0.15; %How de-saturated or transparent to make the patch
if transparent
    faceAlpha=patchSaturation;
    patchColor=col;
    set(gcf,'renderer','openGL')
else
    faceAlpha=1;
    patchColor=col+(1-col)*(1-patchSaturation);
    set(gcf,'renderer','painters')
end

    
%Calculate the y values at which we will place the error bars
uE=errBar(1,:);
lE=errBar(2,:);



%Add the error-bar plot elements
holdStatus=ishold;
if ~holdStatus, hold on,  end


%Make the cordinats for the patch
yP=[lE,fliplr(uE)];
xP=[x,fliplr(x)];

%remove any nans otherwise patch won't work
xP(isnan(yP))=[];
yP(isnan(yP))=[];


H.patch=patch(xP,yP,1,'facecolor',patchColor,...
              'edgecolor','none',...
              'facealpha',faceAlpha);


%Make nice edges around the patch. 
H.edge(1)=plot(x,lE,'-','color',edgeColor);
H.edge(2)=plot(x,uE,'-','color',edgeColor);

%The main line is now covered by the patch object and was plotted first to
%extract the RGB value of the main plot line. I am not aware of an easy way
%to change the order of plot elements on the graph so we'll just remove it
%and put it back (yuk!)
delete(H.mainLine)
H.mainLine=plot(x,y,lineProps{:});


if ~holdStatus, hold off, end

