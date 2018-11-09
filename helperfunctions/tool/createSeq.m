function seq = createSeq(rawFitResultsCam1,rawFitResultsCam2,coordsCam2,tformTotal,beginFrame)
NFrames  =rawFitResultsCam2.Frame(end)+1-beginFrame;
maskFilt2 = true(size(rawFitResultsCam2.Photons,1),1);
maskFilt1 = true(size(rawFitResultsCam1.Photons,1),1);

photons1 = rawFitResultsCam1.Photons(maskFilt1,:);
bg2 = rawFitResultsCam2.Bg(maskFilt2,:);
photons2 = rawFitResultsCam2.Photons(maskFilt2,:);
coords1 = rawFitResultsCam2.Coord(maskFilt2,:);
coords2 = rawFitResultsCam2.Coord(maskFilt2,:);

photons2I = reshape(photons2,[size(photons1,1)/NFrames NFrames]);
bg2I = reshape(bg2,[size(photons2,1)/NFrames NFrames]);

coords1x = reshape(coords1(:,1),[size(photons1,1)/NFrames NFrames]);
coords1y = reshape(coords1(:,2),[size(photons1,1)/NFrames NFrames]);


coords2x = reshape(coords2(:,1),[size(photons1,1)/NFrames NFrames]);
coords2y = reshape(coords2(:,2),[size(photons1,1)/NFrames NFrames]);

esigmax = reshape(rawFitResultsCam2.Sigma(:,1),[size(photons1,1)/NFrames NFrames]);
esigma = reshape(rawFitResultsCam2.Sigma(:,2),[size(photons1,1)/NFrames NFrames]);

coordsCam2x = reshape(coordsCam2(:,1),[size(photons1,1)/NFrames NFrames]);
coordsCam2y = reshape(coordsCam2(:,2),[size(photons1,1)/NFrames NFrames]);


%%

delta = transformPointsInverse(tformTotal,rawFitResultsCam1.Coord)-rawFitResultsCam2.Coord;
deltax = reshape(delta(:,1),[size(delta(:,1),1)/NFrames NFrames]);
deltay = reshape(delta(:,2),[size(delta(:,2),1)/NFrames NFrames]);
PSFSigma=1.39;
seq=[];
seq(:,:,1) = min(esigma,10);
seq(:,:,2) = min(esigmax,10);
seq(:,:,3) = photons2I;
seq(:,:,4) = coords2x-coordsCam2x;
seq(:,:,4)=seq(:,:,4)-repmat(mean(seq(:,:,4),2),[1 size(seq,2)]);
seq(:,:,5) = coords2y-coordsCam2y;
seq(:,:,5) =seq(:,:,5)-repmat(mean(seq(:,:,5),2),[1 size(seq,2)]);

end