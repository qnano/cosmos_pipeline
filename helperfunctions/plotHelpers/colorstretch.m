function cm = colorstretch(value,min_max,cmap)
% COLORSTRETCH    make color map to stretch across given range
% 
% cm = colorstretch(value,min_max)


if nargin < 2 || isempty(min_max)
    min_max = [500 850];
end
if nargin < 3 || isempty(cmap)
    cmap = rgbstretch(0:100,[0 100]);
end

% if sum(value < min_max(1)) || sum(value > min_max(2))
%     error('value must be within range of min_max')
% end

if min_max(1) == min_max(2)
    cm = cmap(round(length(cmap)/2),:);
else
    value0 = round((value-min_max(1))/diff(min_max)*(size(cmap,1)-1))+1;
    value0(value0 < 1) = 1;
    value0(value0 > size(cmap,1)) = size(cmap,1);
    cm = cmap(value0,:);
end



