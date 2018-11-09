function pc=glimpse2avi(gfolder,frames,frmave,bisect,lim,dispscale1,dispscale2,color1,color2,imglbl,FOVlbl1,FOVlbl2,mag,frmrate,avifolder);

%% Notes
%
% 3/3/13
% This is my attempt to understand and modify a routine (macavi_with_overlay_glimpse_v5) that converts
% GLIMPSE image files to movies. I want to create false colored movies
% (e.g. green image channel displays in green) with titles and labels. 
%
%% Input Variables

% gfolder = alphanumeric path to the GLIMPSE images
% frames =  image vector 
% frmave = the number of frames to average
% bisect = row to bisect image to indepedently scale and color FOVs
% lim = [min max] start and end columns to crop FOV 
% dispscale1 = [min max] sets the scale for first FOV
% dispscale2 = [min max] sets the scale for second FOV
% color1 = selects color for first FOV    1=red, 2=green, 3=blue, 4=gray
% color2 = selects color for second FOV
% imglbl =  text label for movie
% FOVlbl1 = text label for first FOV   text color set by color variable
% FOVlbl2 = text label for second FOV
% frmrate== output frame rate for the avi (in frames per second)
% avifolder == folder in which avi will be placed

% load GLIMPSE header file
eval(['load ' [gfolder 'header.mat'] ' -mat']);

% load first image and get frame size
img =glimpse_image(gfolder,vid,1);           
[a b]=size(img);

% create figure
fig=figure(10);
set(10,'DoubleBuffer','on');
set(10,'Name',gfolder);

movieindx=1;
for inum=frames
    for rdindx=inum:inum+frmave-1
        img(:,:,rdindx-inum+1)=glimpse_image(gfolder,vid,rdindx);
    end
    img=imdivide(sum(img,3) ,frmave);
       
% convert average image from double to integers by rounding
img=round(img);

%split images into 2 FOVs and rotate each 90 degrees
FOV1=img(1:bisect,lim(1):lim(2));
FOV1=rot90(FOV1);
FOV2=img(bisect+1:end,lim(1):lim(2));
FOV2=rot90(FOV2);

% color or grayscale maps
color1=uint8(color1);
cmap1 = zeros(dispscale1(2),3);
if color1 ~= 4
    for i=1:max(max(FOV1))
        if i < dispscale1(1)
            cmap1(i,color1)=0;
        elseif i>=dispscale1(1) &&i < dispscale1(2)
            cmap1(i,color1)=(i-dispscale1(1))/(dispscale1(2)-dispscale1(1));
        else
            cmap1(i,color1)=1;
        end    
    end  
elseif color1 == 4
    for i=1:max(max(FOV1))
        if i < dispscale1(1)
            cmap1(i,:)=0;
        elseif i>=dispscale1(1) &&i < dispscale1(2)
            cmap1(i,:)=(i-dispscale1(1))/(dispscale1(2)-dispscale1(1));
        else
            cmap1(i,:)=1;
        end    
    end 
end
FOV1C=ind2rgb(FOV1,cmap1);

color2=uint8(color2);
cmap2 = zeros(dispscale2(2),3); 
if color2 ~= 4
    for i=1:max(max(FOV2))
        if i < dispscale2(1)
            cmap2(i,color2)=0;
        elseif i>=dispscale2(1) &&i < dispscale2(2)
            cmap2(i,color2)=(i-dispscale2(1))/(dispscale2(2)-dispscale2(1));
        else
            cmap2(i,color2)=1;
        end    
    end 
elseif color2 == 4 % grayscale
  for i=1:max(max(FOV2))
        if i < dispscale2(1)
            cmap2(i,:)=0;
        elseif i>=dispscale2(1) &&i < dispscale2(2)
            cmap2(i,:)=(i-dispscale2(1))/(dispscale2(2)-dispscale2(1));
        else
            cmap2(i,:)=1;
        end    
    end 
end
FOV2C=ind2rgb(FOV2,cmap2);

% recombine images
img_recom=horzcat(FOV1C,FOV2C);

%display and label image
imshow(img_recom, 'InitialMagnification', mag)
fi=strfind(gfolder,'\');
filename=gfolder(fi(end-1)+1:fi(end)-1);
filename=strrep(filename,'_','\_');
title([filename ': ' imglbl],'FontSize',14);

% label FOVs
t1=text(20,20,{FOVlbl1},'FontSize',14);
if color1==1
    set(t1,'color','r');
elseif color1==2
    set(t1,'color','g');  
elseif color1==3
    set(t1,'color','b');
elseif color1==4
    set(t1,'color','w');
end  

t2=text(bisect+20,20,{FOVlbl2},'Color','r','FontSize',14);
if color2==1
    set(t2,'color','r');
elseif color2==2
    set(t2,'color','g');  
elseif color1==3
    set(t2,'color','b');
elseif color1==4
    set(t2,'color','w');
end
text(20,500,['frame # ' num2str(inum)],'Color','w','FontSize',14)

% add the current figure to the movie
pc(movieindx)=getframe(gca);                 
movieindx=movieindx+1;
    if movieindx/20==round(movieindx/20)
       disp(['movieindx= ' num2str(movieindx)]);    
    end
end                                                

%sprintf('making an avi file \n')
movie2avi(pc,avifolder,'fps',frmrate,'quality',100,'compression','none');

