function [ mask ] =  filterFits(rawFitResults,paramsFilterFits,fitSigma)
% FILTERFITS    filter rawFitResults property to populate FitResults property

%% Filter fits

mask1=(mean(rawFitResults.CRLB_STD,2) <= paramsFilterFits.MaxCRLBSTD) &...
    (mean(rawFitResults.CRLB_STD,2) >= paramsFilterFits.MinCRLBSTD);
fprintf('Eliminated %d localisations based on localization precision\n',sum(~mask1));

[ ~,~,PFA_adj ]=fdr_bh(reshape(rawFitResults.LL(3,:),[prod(size(rawFitResults.LL(3,:))) 1]),min(0.05,paramsFilterFits.MaxPValue),'dep','no');
mask2=(PFA_adj <= paramsFilterFits.MaxPValue);
fprintf('Eliminated %d localisations based on probabilty of false positive\n',sum(~mask2));

mask3=(rawFitResults.Photons <= paramsFilterFits.MaxPhotons) &...
    (rawFitResults.Photons >= paramsFilterFits.MinPhotons);
fprintf('Eliminated %d localisations based on photons\n',sum(~mask3));

mask4=(rawFitResults.Bg <= paramsFilterFits.MaxBg) &...
    (rawFitResults.Bg >= paramsFilterFits.MinBg);
fprintf('Eliminated %d localisations based on bg\n',sum(~mask4));

cds = [rawFitResults.Coord(:,2) rawFitResults.Coord(:,1)]-[rawFitResults.RoiStart(:,2) rawFitResults.RoiStart(:,1)];
% throw away fits if outside fit box
mask5=(cds(:,1)' < size(rawFitResults.ROIStack,1))'&(0 < cds(:,1)')'...
    &(cds(:,2)' < size(rawFitResults.ROIStack,2))'&(0 < cds(:,2)')';

fprintf('Eliminated %d localisations based on location\n',sum(~mask5));
mask = mask1.*mask2.*mask3.*mask4.*mask5;
if nargin >2 && fitSigma
mask6=(rawFitResults.Sigma(:,1) <= paramsFilterFits.MaxSigmax) &...
    (rawFitResults.Sigma(:,1) >= paramsFilterFits.MinSigmax) &...
    (rawFitResults.Sigma(:,2) <= paramsFilterFits.MaxSigmay) &...
    (rawFitResults.Sigma(:,2) >= paramsFilterFits.MinSigmay);
    fprintf('Eliminated %d localisations based on Sigmaxy\n',sum(~mask6));
mask = mask.*mask6;
end

%Write preFitResults
mask7 = ones(size(mask));
if (paramsFilterFits.MinPixelDist > 0)
    coordsCam1=[rawFitResults.Coord rawFitResults.Frame];
   [idx,dist] = knnsearch(coordsCam1.*repmat([1 1 2*paramsFilterFits.MinPixelDist],[size(coordsCam1,1) 1]),coordsCam1.*repmat([1 1 2*paramsFilterFits.MinPixelDist],[size(coordsCam1,1) 1]),'k',2);     
mask7(dist(:,2) <=paramsFilterFits.MinPixelDist)=0;
mask7(idx(dist(:,2) <=paramsFilterFits.MinPixelDist,2))=0;
end
fprintf('Eliminated %d localisations smaller than minPixelDist\n',size(mask,1)-sum(mask7))
mask=mask.*mask7;
% %%
% %Initialize FitResults property
% FitResults = struct('xCoord',[],'xSigma',[],'yCoord',[],'ySigma',[],'photons',[],'bg',[],...
%     'merit',[],'roiIdxAll',[],'roiIdxFrame',[],'model','modelIdx');
% eliminated=0;
% %Write FitResults
% for ii = 1:max(rawFitResults.Frame)+1
%     
%     %all rawFitResults in frame ii (1-based)
%     frameMask = (rawFitResults.Frame+1)==ii;
%     maskii = frameMask & mask;
%     copymaskii = maskii;
%     %search for localizations within minPixelDist and select the one with
%     %the lowest CRLB
%     if sum(maskii) > 1
%         findmaskii = find(maskii);
%         for jj = findmaskii'
%             %find distances between localization
%             dist = sqrt((rawFitResults.Coord(jj,1) - rawFitResults.Coord(maskii,1)).^2+...
%                 (rawFitResults.Coord(jj,2) - rawFitResults.Coord(maskii,2)).^2);
%             %find distances less than minPixelDist
%             closLocIdx = find(dist < paramsFilterFits.MinPixelDist);
%             if numel(closLocIdx) > 1
%                 eliminated=eliminated+size(closLocIdx,1);
% %                 warning('double localisations smaller than minPixelDist. Elimirnated!!')
% %                 %find min CRLB for set
% %                 [~, minCRLBidx] = min(rawFitResults.CRLB_STD(findmaskii(closLocIdx),1)+...
% %                     rawFitResults.CRLB_STD(findmaskii(closLocIdx),2));
% %                 %alter maskii
% %                 closLocIdx(minCRLBidx) = [];
%                 maskii(findmaskii(closLocIdx)) = 0;
%             end
%         end
%     end
%     mask(copymaskii)=maskii(copymaskii);
%     % gett all fit rsults for frame ii
%     FitResults(ii).PFA = PFA_adj(maskii);
%     FitResults(ii).PH1 = PFA_adj(maskii);
%     
%     
%     FitResults(ii).xCoord = [rawFitResults.Coord(maskii,1) rawFitResults.CRLB_STD(maskii,1)];
%     FitResults(ii).yCoord = [rawFitResults.Coord(maskii,2) rawFitResults.CRLB_STD(maskii,2)];
%     FitResults(ii).photons = [rawFitResults.Photons(maskii) rawFitResults.Photons_STD(maskii)];
%     FitResults(ii).bg = [rawFitResults.Bg(maskii) rawFitResults.Bg_STD(maskii)];
%     ;
% %     FitResults(ii).merit = rawFitResults.LL(maskii);
%     FitResults(ii).roiIdxAll = find(maskii);
%     FitResults(ii).roiIdxFrame = find(maskii(frameMask));
%     FitResults(ii).model = ones(sum(maskii),1);
%     FitResults(ii).modelIdx = ones(sum(maskii),1);
%     if nargin >2 && fitSigma
%         FitResults(ii).xSigma = [rawFitResults.Sigma(maskii,1) rawFitResults.Sigma_STD(maskii,1)];
%         FitResults(ii).ySigma = [rawFitResults.Sigma(maskii,2) rawFitResults.Sigma_STD(maskii,2)];
%     end   
% end
%fprintf('Eliminated %d localisationss smaller than minPixelDist\n',eliminated)

mask=logical(mask);