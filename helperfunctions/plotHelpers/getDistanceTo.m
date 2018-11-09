function [ trackLogA ] = getDistanceTo( maskPeriphery,closedMask, trackLogA )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    

    for l=1:length(trackLogA)
%         clear distanceSM;
        for k=1:length(trackLogA{l}.data)
            if ~isempty(trackLogA{l}.data{k})
                location = trackLogA{l}.data{k}(:,[3 2]);
                clear distanceSM;
                clear logicIdxInMask;
                distanceSM=zeros(size(maskPeriphery,3),size(location,1));
                distanceSMInner=zeros(size(maskPeriphery,3),size(location,1));
                distanceSMouter=zeros(size(maskPeriphery,3),size(location,1));
                for i=1:size(maskPeriphery,3)
                    outterPixel=bdilation(closedMask(:,:,i-1))-closedMask(:,:,i-1);
                    innerPixel=berosion(closedMask(:,:,i-1))-berosion(closedMask(:,:,i-1),2); 
                    img=squeeze(maskPeriphery(:,:,i-1));

                    numberOfSM=size(location,1);
                    blocksize=3;
                    for j=1:blocksize:(ceil(numberOfSM/blocksize)*blocksize);
                         currectBlockSize = size(location(j:min(j-1+blocksize,size(location,1))),2);

                         %Make grid of distance transform
                        [xsg ysg ] = meshgrid(1:size(img,2),1:size(img,1));

                        X=repmat(xsg,[1 1 currectBlockSize]);
                        Y=repmat(ysg,[1 1 currectBlockSize]);

                        %Distance transform of a point
                        R=sqrt((X-repmat(permute(location(j:min(j-1+currectBlockSize,size(location,1)),1),[3 2 1]),[size(img) 1])).^2 ...
                               +(Y-repmat(permute(location(j:min(j-1+currectBlockSize,size(location,1)),2),[3 2 1]),[size(img) 1])).^2 ...
                           );

                        %calucalte distance to nucleus
                        distanceSM(j:j+size(R,3)-1,i) = double(squeeze(min(dip_image(double(repmat(R,[1 1 1 size(img,3)]))),logical(permute(repmat(img,[1 1 size(R,3)]),[2 1 3])),[1 2])))';
                        distanceSMInner(j:j+size(R,3)-1,i) = double(squeeze(min(dip_image(double(repmat(R,[1 1 1 size(innerPixel,3)]))),logical(permute(repmat(innerPixel,[1 1 size(R,3)]),[2 1 3])),[1 2])))';
                        distanceSMouter(j:j+size(R,3)-1,i) = double(squeeze(min(dip_image(double(repmat(R,[1 1 1 size(outterPixel,3)]))),logical(permute(repmat(outterPixel,[1 1 size(R,3)]),[2 1 3])),[1 2])))';
                    end
                end
                
%                 t=0;
%                  A = round(location);
%                 for i=1:size(closedMask,3)
%                     for j= 1:size(A,1)
%                        B = newim(squeeze(closedMask(:,:,i-1)),'double');   
%                        B(A(j,2),A(j,1))= 1;
% %                        disp(A(j,1))
% %                        disp(A(j,2))
%                        idx = (B.*(closedMask(:,:,i-1) & ~maskPeriphery(:,:,i-1)));
%                        logicIdxInMask(j,i) = sum(idx);
%                     end
%                 end
             multiplier = ones(size(distanceSMouter));
             idx = (distanceSMInner <  distanceSMouter   );
                multiplier(idx == 1) = -1;
                trackLogA{l}.distanceSM{k} = distanceSM.*multiplier; 
            end
%             maskPeriphery(A(:,2),A(:,1),0)= 1;
        end
    end
end

