function imh = plotRastergramFromCev(ah,initial_nAOIs,nAOIs,shift,fps,width,cev,color1,color2)


%% expand # of rows of matrix for display
cev_display=zeros(size(cev,1)*width,size(cev,2)); % initialize a cev matrix that has a number of rows = width*cev rows
for j=1:size(cev,1)
    cev_display(width*(j-1)+1:j*width,:)=repmat(cev(j,:),[width,1]); % replicate row j of cev width number of times
end
% for cev_display:
% columns = 0 represent no events
% columns = 1 represent events for interval
% columns = 10 shifts of interval
%% colorize rastergram
% delete(imh)
cev_rgb=repmat(cev_display,[1 1 3]); % create z=3 matrix for rgb display

cev_rgb(:,:,1)=color2(1)*(cev_display==10);  
cev_rgb(:,:,2)= color2(2)*(cev_display==10);             % white
cev_rgb(:,:,3)= color2(3)*(cev_display==10);              % white    
           
for i=0:round(max(cev_display(cev_display<10)))
    cev_rgb(:,:,1)=cev_rgb(:,:,1)+color1(i+1,1).*(cev_display==i);  
    cev_rgb(:,:,2)=cev_rgb(:,:,2)+color1(i+1,2).*(cev_display==i);             % white
    cev_rgb(:,:,3)=cev_rgb(:,:,3)+color1(i+1,3).*(cev_display==i);              % white    
end

if round(max(cev_display(cev_display<10))) > 1
    colorbar('peer',ah)
end

% cev_rgb = cev_rgb./max(cev_rgb);
if isempty(ah)
    ah=axes;
end

imx=((1:size(cev_rgb,2))-shift)/fps;  % shift indices in x such that alignment occurs at x = 0 and convert to seconds
imy=(1:size(cev_rgb,1))/width; % divide by width to get AOI count on y axis
imh=image(imx,imy,cev_rgb,'Parent',ah); % display rastergram 
xlabel('Time (s)'); % label x axis
ylabel('Number of spots');  % label y axis
title(['total # of AOIs = ' num2str(initial_nAOIs)]);
end
