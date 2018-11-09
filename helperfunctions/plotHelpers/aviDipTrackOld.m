function [ h, colors ] = aviDipTrack(tracks,trackVal,sequence, pixelSize,TimeStep,movieFileName, xzoom, yzoom,tzoom)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    % Optional
    % xzoom = []; %make empty for full image
    % yzoom = []; %make empty for full image
    % vidObj.FrameRate = round(1/TimeStep);

    %% example for making an avi
    % pixelSize = 1;
    % TimeStep=1;
    % path = 'C:\Users\smithc5\Dropbox\grunwaldLabPc\summary\detection\'; 
    % movieFileName = fullfile(path,'dipTrackTstMovie.avi');
    % addpath('C:\Users\smithc5\Dropbox\grunwaldLabPc\summary\detection\backup\helperfunctions')
    % fileName{1} = 'LoGforCarlasTH5000.trk';
    % pathName{1} = 'C:\Users\smithc5\Dropbox\grunwaldLabPc\summary\detection\backup\localize\';
    % [ ~,~,tracks ] = importTrackFile(fileName,pathName, 1 )

    %length of scaleLine in microns
    % scaleLineLength = 1; 

    if nargin < 8 || isempty(xzoom) || isempty(yzoom) || isempty(tzoom)
        xzoom = []; %make empty for full image
        yzoom = []; %make empty for full image
        tzoom=size(tracks,2);
    end
    pixelSize = pixelSize/100;
    %set limits for ROI

    if min(trackVal) == max(trackVal)
        optoins = plotTracksV1SetOptions('colorbar',0,'tracks',tracks,'color',[0 0 1]);
    else
        optoins = plotTracksV1SetOptions('colorbar',0,'tracks',tracks,'tracksVal',trackVal);
    end
    % figure;
    % co = hsv(100);
    % plotTracksV1(optoins)
    % hold on
    % lowresolution=true;
    % for i=1:size(closedNucleus,3)
    % %     [r,c,v] = ind2sub(size(squeeze(closedNucleusEdge(:,:,i-1))),find(permute(squeeze(closedNucleusEdge(:,:,i-1)),[1 2]) ==1 ));
    % %     plot(r,c,'or','MarkerSize',3.75)
    % %     hold on
    %         p(i) = patch(isosurface(repmat(logical(permute(squeeze(closedNucleus(:,:,i-1)),[1 2])),[1 1 500]),0.5));
    %         if lowresolution
    %             reducepatch(p(i),0.1)
    %         end
    %         set(p(i),'FaceColor',co(i,:),'EdgeColor','none','FaceAlpha',0.1);
    %         hold on
    % end
    % hold off
    % view([70 70])
    % grid on
    % xlabel('Pixel')
    % ylabel('Pixel')
    % zlabel('Time')


    dipTrackOptions = dipTrackSetOptions;
    dipTrackOptions.plotTracksOptions = optoins;

    dipTrackOptions.im = sequence;
    dipTrackOptions.plotTracksOptions.colorLineByTrackNum = 1;
    % dipTrack(dipTrackOptions)
    h=dipTrack(dipTrackOptions);

    %resize figure
    diptruesize(h,300)

    %zoom into ROI
    pause(0.1)
    ha = findall(h,'type','axes');
    %set xlimits
    if ~isempty(xzoom)
        xlim(ha,xzoom+[-0.5 0.5])
    end
    if ~isempty(yzoom)
        ylim(ha,yzoom+[-0.5 0.5])
    end

    %add appropriate scaleLine
    delete(findall(h,'tag','scaleLine')) %delete current scaleLine
    % scaleLine(ha,scaleLineLength,pixelSize);

    if (nargin > 7 || ~isempty(movieFileName))
        clear mex
        % Prepare VideoWriter object
        vidObj = VideoWriter(movieFileName);
        vidObj.FrameRate = round(TimeStep);
        %open avi file
        open(vidObj);

        for ii = 1:min(size(tracks,2),tzoom)
            %update frame
            dipmapping(h,'slice',ii-1)
            % Write each frame to the file.
            currFrame = getframe(ha);
            data = uint8(extend(currFrame.cdata,[194 194 3]));
            currFrame.cdata = data;
            writeVideo(vidObj,currFrame);
        end

        %close avi file
        close(vidObj);
        close(h);
    end


    minMax = [min(dip_image(trackVal)) max(dip_image(trackVal)) ];
    colors = colorstretch(trackVal(:,1),minMax,optoins.cmap);
end

