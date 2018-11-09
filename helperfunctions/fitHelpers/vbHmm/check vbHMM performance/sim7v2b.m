load('datatestf7.mat')
clear labelest modelHMMest evidenceest modelGMMest indexes
rng('shuffle');

N=64; 
BStrap = [0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.4 0.6 0.8 1];

X(repmat(permute(states == 0,[3 1 2]),[2 1 1])) = NaN;
order = 1:4;
if size(X,3) > 1 
    dataIdx = ~isnan(X);
    X2 = reshape(X(dataIdx),[size(X,1) sum(sum(sum(dataIdx)))./size(X,1)]);
end
prior.m = repmat(mean(mean(X2,2),3),[1  max(order)]);
% prior.m = [200 200 200 200; 100 300 600 900];
prior.alpha = 1;
prior.kappa = 1;
prior.v = 1.*(size(X,1)+1); 
prior.alphaA = 1/2;
prior.M = 1;

%%
for j=1:length(BStrap)
    disp(['j:' num2str(j) 'of ' num2str(length(BStrap))])
    parfor i=1:N
        idx = randperm(size(X,2),round(size(X,2)*BStrap(j)));
        while sum(sum(sum(isnan(squeeze(X(:,idx,:)))))) == numel(X(:,idx,:))
            idx = randperm(size(X,2),round(size(X,2)*BStrap(j)));
        end
        indexes{i} = idx;
        [labelest{i}, modelHMMest{i}, evidenceest{i}, modelGMMest{i}] = bayesHMM(X(:,idx,:), 4,3,prior,0)
        
       disp(['i:' num2str(i) 'of ' num2str(N)])
    end
    save(['HMMloc7dfbl' num2str(j)],'-v7.3')
    clear labelest modelHMMest evidenceest modelGMMest indexes
end

exit