function [label, modelHMM, evidence, modelGMM] = bayesHMMv2(data, m,n, prior,inovation)
        evidence=[];
        priorCollection=[];
        if isscalar(m)
            order = 1:m;
        else
            order =m;
        end
        if nargin < 5 || isempty(inovation)
            inovation = 0.01; %% regularisation to prevent states with very low occupancy.
        end
        
        if nargin > 3
            if isscalar(prior.m) || isvector(prior.m)
                priorCollection.m = repmat(prior.m,[1 max(order)]);
            else
                priorCollection.m=prior.m;
            end
            if isscalar(prior.alpha)
                priorCollection.alpha = repmat(prior.alpha,[1 max(order)]);
            else
                priorCollection.alpha=prior.alpha;
            end
            if isscalar(prior.kappa)
                priorCollection.kappa = repmat(prior.kappa,[1 max(order)]);
            else
                priorCollection.kappa=prior.kappa;
            end
            if isscalar(prior.v)
                priorCollection.v = repmat(prior.v,[1 max(order)]);
            else
                priorCollection.v=prior.v;            
            end
            if isscalar(prior.alphaA) 
                for k=order
                    priorCollection.alphaA{k} = prior.alphaA.*ones(k);
                end
            else
                priorCollection.alphaA{1}=prior.alphaA;    
            end

            if isscalar(prior.M)    
                for k=order
                    priorCollection.M{k} = prior.M.*repmat(eye(size(data,1)),[1 1 k]);   % M = inv(W) 
                end
            else
                priorCollection.M{1}=prior.M;    
            end
        else
               if size(data,3) > 1 
                    dataIdx = ~isnan(data);
                    X = reshape(data(dataIdx),[size(data,1) sum(sum(sum(dataIdx)))./size(data,1)]);
                end
                priorCollection.m = repmat(mean(mean(X,2),3),[1  max(order)]);
                priorCollection.alpha = ones(1,max(order));
                priorCollection.kappa = ones(1,max(order));
                priorCollection.v = ones(1,max(order)).*(size(data,1)+1); 
                for k=order
                    priorCollection.alphaA{k} = ones(k)/2;
                    priorCollection.M{k} = repmat(eye(size(data,1)),[1 1 k]);
                end
        end

        
        for k=order
            for j=1:n
                modelPrior2.alpha = priorCollection.alpha(k);
                modelPrior2.alphaA = priorCollection.alphaA{k};
                modelPrior2.kappa = priorCollection.kappa(k);
                modelPrior2.m = priorCollection.m(:,1:k);
                modelPrior2.v = priorCollection.v(k);
                modelPrior2.M = priorCollection.M{k}
                [y1{k,j}, modelTemp1{k,j}, L] = mixGaussVb(data, k, modelPrior2);
                evidence(k,j) = L(end).*all(round(sum(modelTemp1{k,j}.R,1)./size(modelTemp1{k,j}.R,1)*100)/100 > inovation);      
                
                %*sum(modelTemp1{k,j}.alpha./sum(modelTemp1{k,j}.alpha) > inovation)/length(modelTemp1{k,j}.alpha);
            end
        end

        [~,kopt] = max(mean(evidence,2));
        kopt = kopt(1);
        disp(['Order ' num2str(kopt)])
        plot(mean(evidence,2))
                
        zeroStates = 1;
        while zeroStates
            [~,jopt] = max(evidence(kopt,:)); 
            modelGMM = modelTemp1{kopt(1),jopt(1)};

            modelPrior.alpha = priorCollection.alpha(kopt(1));
            modelPrior.alphaA = priorCollection.alphaA{kopt(1)};
            modelPrior.kappa = priorCollection.kappa(kopt(1));
            modelPrior.m = priorCollection.m(:,1:kopt(1));
            modelPrior.v = priorCollection.v(kopt(1));
            modelPrior.M = priorCollection.M{kopt(1)};

            [label, modelHMM,L] = multiGaussVbv2(data, kopt, modelPrior);
            plot(L)
            zeroStates  = ~all(round(sum(modelHMM.R,1)./size(modelHMM.R,1)*100)/100 > inovation);
            if kopt == 1
                break;
            elseif zeroStates
                kopt=kopt-1;
            end
        end
        if exist('prior')
            modelHMM.modelPrior=prior;
        end
end
