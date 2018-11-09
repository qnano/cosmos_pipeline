function [img,  stackLength] = testShowSlice(params,handle,newSlice)

if nargin == 1
    newSlice=1;
    handle=[];
end


if ndims(params.rgbImage) <= 3
    stackLength = max(1,size(params.rgbImage,3));
    img = repmat(params.rgbImage(:,:,newSlice),[1 1 3]);
elseif ndims(params.rgbImage) == 4 && size(params.rgbImage,3) == 3
    stackLength = max(1,size(params.rgbImage,4));
    img = params.rgbImage(:,:,:,newSlice);
end 

if ~isempty(handle) & ishandle(handle)
    if ndims(params.rgbImage) <= 3
        handle.CData = repmat(params.rgbImage(:,:,newSlice),[1 1 3]);
    elseif ndims(params.rgbImage) == 4 && size(params.rgbImage,3) == 3
        handle.CData = params.rgbImage(:,:,:,newSlice);
    end 
end

end
