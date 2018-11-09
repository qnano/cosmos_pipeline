clear all
close all

addpath(genpath('../../../helperfunctions'))

load('testmodel','modelGMM','modelHMM')
N=0;
%%
NumTraj=400;
Timesteps=1500;
buffer = round(Timesteps*0.1);
d=2;
intergrationTime = 5;

% ton = 10;
% toff = 5;
ton = 25;
toff = 5;

% TransMatrix1 = [1-(ton*intergrationTime)^-1 (ton*intergrationTime)^-1 (ton*intergrationTime)^-1; (ton*intergrationTime)^-1 1-(ton*intergrationTime)^-1 (ton*intergrationTime)^-1;    (toff*intergrationTime)^-1 (toff*intergrationTime)^-1 1-(toff*intergrationTime)^-1];
TransMatrix1 = [1-(ton*intergrationTime)^-1 (ton*intergrationTime)^-1; (toff*intergrationTime)^-1 1-(toff*intergrationTime)^-1;];

TransMatrix1 = mk_stochastic(TransMatrix1);

[states] = genHMMseq(NumTraj, Timesteps, TransMatrix1,ones(NumTraj,1)');
esttrans(states(:,:))
TransMatrix1
BStrap = 0.1:0.1:1;
Nsamps=32
for n=1:Nsamps
    for j=1:length(BStrap)
        idx = randperm(size(states,1),round(size(states,1)*BStrap(j)));
        [tr(:,:,j,n)] = esttrans(states(idx,:));
    end
end
%%
figure
pm = mean(tr,4);
ps = std(tr,0,4);
m = toff; 
errorbar(squeeze(pm(1,2,:)),squeeze(ps(1,1,:)))
hline(TransMatrix1(1,2))
hold on
m = ton;
errorbar(squeeze(pm(2,1,:)),squeeze(ps(2,2,:)))
hline(TransMatrix1(2,1))


k=max(states(:));
for i=1:k
    numObs(i) = sum(states(:)==i);
end

n = numel(states);

prior.alpha = modelHMM.alpha;  % hyperparameter of Dirichlet prior
prior.alphaA =modelHMM.alphaA;
for i=1:k
    prior.W0(:,:,i) = 1e3.*eye(2,2);
end
prior.v = modelHMM.v;  % hyperparameter of inverse Wishart prior of covariances
prior.m = modelHMM.m;  % hyperparameter of Guassian prior of means
prior.kappa = modelHMM.kappa;
w = dirichletRnd(prior.alpha',ones(1,k)/k);
z = states;

mu = zeros(d,k);
Sigma = zeros(d,d,k);
X = zeros(d,NumTraj,Timesteps);
prior.kappa = 10*ones(k,1); 

prior.m = [200 200 200 200; 100 500 1000 1500];
sw = load('dataf6')
mu = sw.mu;
Sigma = sw.Sigma;

for i = 1:k
    idx = z==i;
%     Sigma(:,:,i) = iwishrnd(prior.W0(:,:,i),d+1).*eye(2,2); % invpd(wishrnd(W0,v0));
%     mu(:,i) = gaussRnd(prior.m(:,i),prior.kappa(i)*Sigma(:,:,i));
    X(:,idx) = gaussRnd(mu(:,i),Sigma(:,:,i),sum(idx(:)));
end
model.mu = mu;
model.Sigma = Sigma;
model.w = w;
figure
for i = 1:k
    idx = z==i;
    scatter(X(1,idx),X(2,idx))
    hold on;
end


%%
% 
% load('datatestf')

% 
% sw = load('datatestf2')
% figure
% for i = 1:k
%     idx = z==i;
%     scatter(X(1,idx),X(2,idx))
%     hold on;
% end
% 
% load('datatestf3')
figure
for i = 1:k
    idx = z==i;
    scatter(X(1,idx),X(2,idx))
    hold on;
end
% 
% plot(squeeze(X(2,10,:))); hold on; plot(0.5.*max(X(2,10,:))*idx(10,:))
% 

scatter(X(1,:),X(2,:))
%%
load('datatestf7')
prior.m = [200.*ones(5,1) 200.*cumsum(ones(5,1));]'; 
prior.m = [200 200 200 200; 200 400+N 600+N 800+N];
X(repmat(permute(states == 0,[3 1 2]),[2 1 1])) = NaN;
order = 1:4;
if size(X,3) > 1 
    dataIdx = ~isnan(X);
    X2 = reshape(X(dataIdx),[size(X,1) sum(sum(sum(dataIdx)))./size(X,1)]);
end
% prior.m = repmat(mean(mean(X2,2),3),[1  max(order)]);
prior.m = [200 200 200 200; 100 500 1000 1500];
prior.alpha = 1;
prior.kappa = 1;
prior.v = 1.*(size(X,1)+1); 
prior.alphaA = 1/2;
prior.M = 1;


[labelest, modelHMMest, evidenceest, modelGMMest] = bayesHMM(X(:,1:100,:), max(order),3,prior,0.05)

data =X(:,1:100,:);    
data = permute(data,[1 3 2]);

dataIdx = ~isnan(data);
a = reshape(data(dataIdx),[2 47718]);

figure
scatter(a(1,labelest==1),a(2,labelest==1))
hold on
scatter(a(1,labelest==2),a(2,labelest==2))

% X1=X(:,1:50,:);
% dataIdx = ~isnan(X1);
% X2 = reshape(X1(dataIdx),[size(X,1) sum(sum(sum(dataIdx)))./size(X,1)]);
% 
% temp = permute(X(:,1:50,:),[1 3 2]);
% dataIdx = ~isnan(temp);
% temp2 = temp(dataIdx);
% temp3 = reshape(temp2,[2 24177])
% for i = 1:max(labelest)
%     idx = labelest ==i;
%     scatter(temp(1,idx),temp(2,idx))
%     hold on;
% end


esttrans(states(1:50,:))
a = zeros(size(squeeze(dataIdx(1,:,:))));
a(states(:,1:50,:) > 0) =labelest;
a(dataIdx(1,:,:)) = labelest; 
save('dataf7')
% esttrans(reshape(labelest,50,1500))