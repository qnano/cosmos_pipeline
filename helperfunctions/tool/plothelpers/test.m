close all
clear all



h = figure;
% uip = uipanel('parent',h,'position',[0 0 0.5 0.5])
uip1 = subplot(1,2,1)
b = imtool3D(uip1,rand(100,100,3,100))
b.haxis 

uip2 = subplot(1,2,2)
a = imtool3D(uip2,rand(100,100,3,100))
a.haxis 
% params.rgbImage = rand(100,100,100);
% imtool3D(uip,@testShowSlice,params)
