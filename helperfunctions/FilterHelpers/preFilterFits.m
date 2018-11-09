function [ mask ] =  preFilterFits(coordsCam1,detParCam1,paramsPreFilterFits,keepone)
% FILTERFITS    filter rawFitResults property to populate FitResults property
%
% Created by Carlas Smith September 2015 (UMASS/TU-Delft)

%% Filter fits
if ~isempty(detParCam1)
    mask1=(detParCam1.circularity <= paramsPreFilterFits.circularityMax) &...
        (detParCam1.circularity >= paramsPreFilterFits.circularityMin);
    fprintf('Eliminated %d localisations based on circularity\n',sum(~mask1));

    mask2=  ((detParCam1.PH1 <= paramsPreFilterFits.PH1Max) &...
        (detParCam1.PH1 >= paramsPreFilterFits.PH1Min))';
    fprintf('Eliminated %d localisations based on detection probability\n',sum(~mask2));

    mask3=(detParCam1.clusterSize <= paramsPreFilterFits.clusterSizeMax) &...
        (detParCam1.clusterSize >= paramsPreFilterFits.clusterSizeMin);
    fprintf('Eliminated %d localisations based on cluster size\n',sum(~mask3));

    mask=mask1.*mask2.*mask3;
else
    mask = ones(size(coordsCam1,1),1);
end

%%
eliminated=0;
%Write preFitResults
mask4 = ones(size(mask));
if (paramsPreFilterFits.minPixelDist > 0)
 [idx,dist] = knnsearch(coordsCam1.*repmat([1 1 2*paramsPreFilterFits.minPixelDist],[size(coordsCam1,1) 1]),coordsCam1.*repmat([1 1 2*paramsPreFilterFits.minPixelDist],[size(coordsCam1,1) 1]),'k',2);   
%     for ii = 1:max(coordsCam1(:,3))+1
% 
%         %all rawFitResults in frame ii (1-based)
%         frameMask = (coordsCam1(:,3)+1)==ii;
%         maskii = frameMask & mask;
%         copymaskii = maskii;
%         %search for localizations within minPixelDist and select the one with
%         %the lowest CRLB
%         if sum(maskii) > 1
%             findmaskii = find(maskii);
%             for jj = findmaskii'
%                 %find distances between localization
%                 dist = sqrt((coordsCam1(jj,1) - coordsCam1(maskii,1)).^2+...
%                     (coordsCam1(jj,2) - coordsCam1(maskii,2)).^2);
%                 %find distances less than minPixelDist
%                 closLocIdx = find(dist < paramsPreFilterFits.minPixelDist);
%                 if numel(closLocIdx) > 1
%     %                 warning('double localisations smaller than minPixelDist. Eliminated!!')
%                     eliminated=eliminated+size(closLocIdx,1);
%                     %alter maskii
%                     maskii(findmaskii(closLocIdx)) = 0;
%                 end
%             end
%         end
%         mask(copymaskii)=maskii(copymaskii);
%     end
mask4(dist(:,2) <=paramsPreFilterFits.minPixelDist)=0;
mask4(idx(dist(:,2) <=paramsPreFilterFits.minPixelDist,2))=0;

end
fprintf('Eliminated %d localisations smaller than minPixelDist\n',size(mask,1)-sum(mask4))

mask=mask&mask4;