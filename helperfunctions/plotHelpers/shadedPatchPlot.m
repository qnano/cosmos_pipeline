function H=shadedPatchPlot(x,Q2,Q,lineProps,patchSaturation)
%shadedPatchPlot Creates DStorm Simlation
% SYNOPSIS:
%    H=shadedPatchPlot(x,Q2,Q,lineProps,patchSaturation)
% 
% PARAMETERS:
%   x: x values
%   Q2: main solid line (mean)
%   Q: cell with patch values 2 vector with upper and lower value
%   lineProps
%   patchSaturation
% 
% OUTPUTS:
%   H Figure handle
% EXAMPLE
% 
% h2 = shadedPatchPlot(t(1:end-1), Q2(2,:), ...
%     {[Q1(2,:); Q3(2,:)],[NO1(2,:); NO2(2,:)]},{'r','Linewidth',2});
% 
%   SEE ALSO:
%       plot
holdStatus = ishold;

Q2=Q2(:)';

if isempty(x)
    x=1:length(Q2);
else
    x=x(:)';
end

if length(x) ~= length(Q2)
    error('inputs x and y are not of equal lengths')
end

defaultProps={'-k'};
if nargin<4 || isempty(lineProps)
    lineProps=defaultProps; 
end
if ~iscell(lineProps)
    lineProps={lineProps}; 
end

% Main line
H.mainLine=plot(x,Q2,lineProps{:});
col=get(H.mainLine,'color');

%draw paches + edges
if nargin < 5
    if length(Q) == 1
        patchSaturation = 0.2;
    elseif length(Q) == 2
        patchSaturation = [0.2, 0.2*0.4];
    else
        patchSaturation= logspace(log10(0.2),log10(0.2^2),length(Q));
    end
end

for ll=1:length(Q)
    Q13=Q{ll};
    if length(x) ~= length(Q13)
        error('inputs x and y must have the same length as Q13')
    end

    [H.patch{ll}, uE, lE] = drawpatch(x,col,patchSaturation(ll),Q13);
    edgeColor = col+(1-col)*0.25;
    H.edge{ll} = drawEdge(x,lE,uE,edgeColor);
end

delete(H.mainLine);
H.mainLine=plot(x,Q2,lineProps{:});

if ~holdStatus, hold off, end
end

function edge = drawEdge(x,lE,uE,edgeColor)
    edge(1)=plot(x,lE,'-','color',edgeColor);
    edge(2)=plot(x,uE,'-','color',edgeColor);
end

function [Hpatch, upE, loE]= drawpatch(x,col,patchSaturation,twoVector)
    faceAlpha=patchSaturation;
    patchColor=col;
    set(gcf,'renderer','openGL')
    upE=twoVector(1,:);
    loE=twoVector(2,:);
    if ~ishold, hold on,  end
    yP=[loE,fliplr(upE)];
    xP=[x,fliplr(x)];
    %remove any nans otherwise patch won't work
    xP(isnan(yP))=[];
    yP(isnan(yP))=[];

    Hpatch=patch(xP,yP,1,'facecolor',patchColor,...
                  'edgecolor','none',...
                  'facealpha',faceAlpha);    
end    
