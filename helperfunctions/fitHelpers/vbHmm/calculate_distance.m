function KLdist = calculate_distance(SigmaEst,muEst,sigmaTrue,muTrue)
for l=1:size(SigmaEst,3) %max state 
    for i=1:size(sigmaTrue,3) %num true states
        mu0 = muTrue(:,i);
        sigma0 = sigmaTrue(:,:,i);

        sigma1 = SigmaEst(:,:,l);
        mu1 = muEst(:,l);
        ndims = length(mu0);
        
        if all(mu1 ~= 0)
            tmp = inv(sigma0)*sigma1;
            dKL(i,l) = 1/8*(mu1-mu0)'*inv(sigma0)*(mu1-mu0);%0.5*(trace(tmp)+(mu1-mu0)'*inv(sigma0)*(mu1-mu0)-ndims-log(det(tmp))); %1/8*(mu1-mu0)'*inv(sigma0)*(mu1-mu0)+1/2*log(det(sigma0+sigma1)/sqrt(det(sigma0)*det(sigma0))); 
        else
            dKL(i,l) = NaN;
        end
    end
end
dKLsort = sort(dKL,1);
KLdist = dKLsort(1,:);
end