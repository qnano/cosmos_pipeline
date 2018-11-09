function [Bp,B2,z] = multiMixGaussVbPred(model,X)

    modelTemp.alpha =  reshape(model.alphaPhi,[size(model.alphaPhi,1)*size(model.alphaPhi,2) size(model.alphaPhi,3)])';
    modelTemp.U = reshape(model.U,[size(model.U,1) size(model.U,2) size(model.U,3)*size(model.U,4)]); 
    modelTemp.v = reshape(model.v,[size(model.v,1)*size(model.v,2) size(model.v,3)])';  
    modelTemp.m = reshape(model.m,[size(model.m,1) size(model.m,2)*size(model.m,3)]); 
    modelTemp.logW = reshape(model.logW,[size(model.logW,1)*size(model.logW,2) size(model.logW,3)])';  
    modelTemp.kappa = reshape(model.kappa,[size(model.kappa,1)*size(model.kappa,2) size(model.kappa,3)]); 

    [~, B] = mixGaussVbPred(modelTemp, X);
    B2 = reshape(B,[size(B,1) size(model.m,2) size(model.m,3)]);
    if isfield(model,'mixMat')
        mix = mk_stochastic(model.mixMat);
    else
       mix=ones(size(model.m,2),size(model.m,3));
       mix = mk_stochastic(mix);
    end
    for i=1:size(model.m,2) 
        for j=1:size(model.m,3)
          Bp(:,i,j) = mix(i,j).*B2(:,i,j);
        end
    end
    
    Bp = sum(Bp,3);
    [~,z(:)] = max(Bp,[],2);
end