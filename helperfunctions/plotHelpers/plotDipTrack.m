function [ h, colors ] = plotDipTrack(tracks,tracksVal,closedNucleus,lowresolution)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if min(tracksVal) == max(tracksVal)
        optoins = plotTracksV1SetOptions('colorbar',0,'tracks',tracks,'color',[0 0 1]);
    else
        optoins = plotTracksV1SetOptions('colorbar',0,'tracks',tracks,'tracksVal',tracksVal);
    end
    co = hsv(10);
    h= plotTracksV1(optoins)
    hold on
    for i=1:size(closedNucleus,3)
%             tracks = dip_image(tracks);
%             tracksx = tracks(:,:,0);
%             tracksy = tracks(:,:,1);
%             indx = floor(min(tracksx(tracksx  >0))):ceil(max(tracks(:,:,0)));
%             indy =floor(min(tracksy(tracksy  >0))):ceil(max(tracks(:,:,1)));
%             [X, Y] = meshgrid(indx,indy);
%             a= sub2ind([size(closedNucleus,1) size(closedNucleus,2)],reshape(X,[1 prod(size(X))]),reshape(Y, [1 prod(size(Y))]));
%             img = zeros(size(closedNucleus));
%             img(reshape(X,[1 prod(size(X))]),reshape(Y, [1 prod(size(Y))]),1) = 1;
            p(i) = patch(isosurface(repmat(logical(permute(squeeze(closedNucleus(:,:,i-1)),[1 2])),[1 1 size(tracks,2)]),0.5));
            if lowresolution
                reducepatch(p(i),0.1)
            end
            set(p(i),'FaceColor',co(i,:),'EdgeColor','none','FaceAlpha',0.1);
            hold on
    end
    hold off
    view([13 80])
    grid on
    xlabel('Pixel')
    ylabel('Pixel')
    zlabel('Time')

    minMax = [min(dip_image(tracksVal)) max(dip_image(tracksVal)) ];
    colors = colorstretch(tracksVal(:,1),minMax,optoins.cmap);
end

