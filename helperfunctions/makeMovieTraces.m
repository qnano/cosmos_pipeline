function makeMovieTraces(traces2_allfr,seq,rawFitResultsCam1,rawFitResultsCam2,coordsCam2,movieFileName2,Sample,TH)

if nargin < 7
    Sample = 1;
end

if nargin < 8
    TH = [];
end



NFrames  =rawFitResultsCam2.Frame(end)+1;
maskFilt2 = true(size(rawFitResultsCam2.Photons,1),1);
maskFilt1 = true(size(rawFitResultsCam1.Photons,1),1);

photons1 = rawFitResultsCam1.Photons(maskFilt1,:);
frames2 = rawFitResultsCam2.Frame(maskFilt2,:);
frames1 = rawFitResultsCam1.Frame(maskFilt1,:);
bg2 = rawFitResultsCam2.Bg(maskFilt2,:);
photons2 = rawFitResultsCam2.Photons(maskFilt2,:);
% coords2 = rawFitResultsCam2.Coord(maskFilt2,:);

photons2I = reshape(photons2,[size(photons1,1)/NFrames NFrames]);
bg2I = reshape(bg2,[size(photons2,1)/NFrames NFrames]);

photons1I = reshape(photons1,[size(photons1,1)/NFrames NFrames]);
% coords2x = reshape(coords2(:,1),[size(photons1,1)/NFrames NFrames]);
% coords2y = reshape(coords2(:,2),[size(photons1,1)/NFrames NFrames]);
% 
% esigmax = reshape(rawFitResultsCam2.Sigma(:,1),[size(photons1,1)/NFrames NFrames]);
% esigma = reshape(rawFitResultsCam2.Sigma(:,2),[size(photons1,1)/NFrames NFrames]);
% 
% coordsCam2x = reshape(coordsCam2(:,1),[size(photons1,1)/NFrames NFrames]);
% coordsCam2y = reshape(coordsCam2(:,2),[size(photons1,1)/NFrames NFrames]);
% 

NFrame = rawFitResultsCam1.Frame(end)+1;

visSpots = randperm(size(photons1,1)/NFrame,min(size(photons1,1)/NFrame,25));


% traces2_allfr=~(path-1);
h = figure('position', [100, 100, 800, 700])

%traces2_allfr = sqrt((coords2x-coordsCam2x ).^2+ (coords2y-coordsCam2y ).^2) < 1.5 & photons2I > 10*bg2I;
% scaleLineLength=1;
diptruesize(h,200)
% pixelSize=0.1;


% Prepare VideoWriter object
TimeStep=15;
vidObj = VideoWriter(movieFileName2);
vidObj.FrameRate = round(TimeStep);

%open avi file
open(vidObj);
PSFSigma=1.4;
% visSpots  = randi(size(photons1,1)/NFrame,25,1);
for jj=1:length(visSpots)
Nspot = visSpots(jj);
    for ii = 0:Sample:NFrame-1

        idx = find(frames2==ii);
        A(:,:,ii+1) = rawFitResultsCam2.ROIStack(:,:,idx(Nspot));
        idx1 = find(frames1==ii);
        A1(:,:,ii+1) = rawFitResultsCam1.ROIStack(:,:,idx1(Nspot));
        xxyy(ii+1,:) = rawFitResultsCam2.Coord(idx(Nspot),:)-rawFitResultsCam2.RoiStart(idx(Nspot),:);
        xxyytarget(ii+1,:) = coordsCam2(idx(Nspot),1:2)-rawFitResultsCam2.RoiStart(idx(Nspot),:);

        subplot(2,1,1)

        plot(photons2I(Nspot,:),'-xr')
        hold on
        plot(photons1I(Nspot,:),'-g')
        plot(bg2I(Nspot,:)*(PSFSigma*2)^2,'Color',[.7 .5 0])
        if ~isempty(TH)
            hline(TH, '-k')
%             plot(100*photons2I(Nspot,:)./sqrt(photons2I(Nspot,:).^2+(bg2I(Nspot,:)*(PSFSigma*2)^2).^2),'xk')
        end
        plot(traces2_allfr(Nspot,:)*max(max(photons2I(Nspot,:)),max(photons1I(Nspot,:))),'-b')
        ylabel('I [# Photons]')
        xlabel('Time [s]')
        axis tight
        legend('Complex','Target','Background','Binary','Location','southoutside','Orientation','horizontal')
        
        xlim([1, size(photons2I,2)]);
        if size(seq,3) > 3
            b = sqrt(seq(jj,:,4).^2+seq(jj,:,5).^2);
        end
        I=seq(jj,:,3);
        
        if size(seq,3) > 3
            title(sprintf('spot ID = %d\n mean(\\Delta X) = %0.2g, std(\\Delta X) = %0.2g, max(\\Delta X) = %0.2g, min(I) = %0.2g',Nspot,mean(b(~(traces2_allfr(jj,:)-1))),std(b(~(traces2_allfr(jj,:)-1))),max(b(~(traces2_allfr(jj,:)-1))),min(I(~(traces2_allfr(jj,:)-1)))));
        end
        vline(ii,'-k','t')
        hold off
        hs = subplot(2,1,2);

        axes1Position = get(gca, 'Position');
        delete(hs);
        
        % Position the logo in the upper right.
        x1=0;
        y1=0.01;
        hAxis2 = axes('Position', [x1 y1 0.5 0.5]);
        axis off; % Turn off tick marks, etc.
        gimr = mat2gray(A(:,:,ii+1)'); %flip(flip(A(:,:,ii+1),1),2));
        imshow(gimr );
        hold on
        plot(xxyy(ii+1,1)+1,xxyy(ii+1,2)+1,'h','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',8)
        plot(xxyytarget(ii+1,1)+1,xxyytarget(ii+1,2)+1,'h','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',8)
        title('Complex')
        hAxis2 = axes('Position', [x1+0.5 y1 0.5 0.5]);
        
        gimg = mat2gray(flip(flip(A1(:,:,ii+1),1),2)); 
        gimg = uint8(gimg * 256);
        gimr = uint8(gimr * 256);
        filler = zeros(size(gimr),'uint8');
        rgbImage = cat(3,gimr,gimg,filler);


        imshow(rgbImage);
        hold on
        plot(-10,-10,'h','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8)
        plot(xxyy(ii+1,1)+1,xxyy(ii+1,2)+1,'h','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',8)
        plot(xxyytarget(ii+1,1)+1,xxyytarget(ii+1,2)+1,'h','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',8)
                
        legend('Position')
        title('Target & Complex')
        currFrame = getframe(h);
        data = uint8(extend(currFrame.cdata,[194 194 3]));
        currFrame.cdata = data;
        writeVideo(vidObj,currFrame);

    end
end

%close avi file
close(vidObj);
close(h);
