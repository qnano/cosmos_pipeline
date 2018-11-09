function [hl ht] = scaleLine(ha,len,pxlsize,loc,color,width,fontsize)
% SCALELINE     make scalebar for axes
% 
% [hl ht] = scaleLine(ha,len,pxlsize,loc,color,width,fontsize)
% 
% INPUTS
%   ha - axes handle. Default gca.
%   len - desired length of scale bar um. Default 5.
%   pxlsize - size of pixels in um. Default 0.267.
%   loc - desired location [x y]. x specifies the right edge of the 
%         scalebar (Default sum(max(get(ha,'xlim'))*[1 -.1])). y specifies 
%         the y position of the line sum(max(get(gca,'ylim'))*[1 -.1]).
%   color - string or vector of length 3 indicating the desired color of
%           the scalebar. Default [1 1 1].
%   width - desired line width. Default 3.
%   fontsize - desired fontsize. Default 16.
% OUPUTS
%   hl - line handle
%   ht - text handle



%% Default values
if nargin <1 || isempty(ha)
    ha = gca;
end
if nargin <2 || isempty(len)
    len = 5;
end
if nargin <3 || isempty(pxlsize)
    pxlsize = 0.267;
end
xlim = get(ha,'xlim');
ylim = get(ha,'ylim');
ydist = (ylim(2)-ylim(1))/200;
xdist = (xlim(2)-xlim(1))/100;
if nargin <4 || isempty(loc)
    loc = [max(xlim) max(ylim)];
end
if nargin <5 || isempty(color)
    color = [1 1 1];
end
if nargin <6 || isempty(width)
    width = 3;
end
if nargin <7 || isempty(fontsize)
    fontsize = 11;
end


%% Body
if loc(1) + len/pxlsize > max(xlim)
    loc(1) =   max(xlim)-xdist;
end
ht = text(loc(1)-len/pxlsize/2,loc(2)+ydist,[num2str(len) ' \it\mum'],...
    'verticalalignment','top','horizontalalignment','center','fontweight',...
    'bold','fontsize',fontsize,'color',color,'parent',ha,'tag','scaleLine');
a = get(ht,'extent');
if a(1)+a(3) > max(xlim)
    %     loc(1) =  sum(max(get(ha,'xlim'))*[1 -.03]);
    loc(1) = max(xlim)-a(3)/2;
    delete(ht)
%     ht = text(loc(1)-len/pxlsize/2,loc(2)+ydist,[num2str(len) ' \it\mum'],'verticalalignment','top','horizontalalignment','center','fontweight','bold','fontsize',fontsize,'color',color);
    ht = text(loc(1)-len/pxlsize/2,loc(2)+ydist,[num2str(len) ' \it\mum'],...
        'verticalalignment','top','horizontalalignment','center',...
        'fontweight','bold','fontsize',fontsize,'color',color,'parent',ha,'tag','scaleLine');
end
if a(2)+a(4) > max(ylim)
    %     loc(2) =  sum(max(get(ha,'ylim'))*[1 -.11]);
    loc(2) =   max(ylim)-a(4)-ydist;
    delete(ht)
    ht = text(loc(1)-len/pxlsize/2,loc(2)+ydist,[num2str(len) ' \it\mum'],...
        'verticalalignment','top','horizontalalignment','center',...
        'fontweight','bold','fontsize',fontsize,'color',color,'parent',ha,'tag','scaleLine');
end  
hl = line([loc(1)-len/pxlsize loc(1)],[loc(2) loc(2)],'color',color,...
    'linewidth',width,'parent',ha,'tag','scaleLine');
set(get(ha,'parent'),'inverthardcopy','off')
