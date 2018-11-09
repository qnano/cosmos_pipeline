clear all
close all

addpath(genpath('../../../helperfunctions'))

load('testmodel','modelGMM','modelHMM')
% load('traces','traces2_allfr', 'ttb')
N=0;
%%
NumTraj=400;
Timesteps=1500;
buffer = round(Timesteps*0.1);
d=2;
intergrationTime = 5;

% ton = 10;
% toff = 5;
ton = 16;
toff = 8;
TransMatrix1 = [1-(ton*intergrationTime)^-1 (ton*intergrationTime)^-1;  (toff*intergrationTime)^-1 1-(toff*intergrationTime)^-1];
[statesBind] = genHMMseq(NumTraj, Timesteps, TransMatrix1,ones(NumTraj,1)');

BStrap = 0.1:0.1:1;
Nsamps=32
for n=1:Nsamps
    for j=1:length(BStrap)
        idx = randperm(size(statesBind,1),round(size(statesBind,1)*BStrap(j)));
        [tr(:,:,j,n)] = esttrans(statesBind(idx,:));
    end
end

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


%%
align = false;
delete = true;
sortType = 3;
sortState = 1;
colorArray =  [1 1 1; 0 0 1; 1 0 0; 1 0 0; 0 1 0; 1 0.4 0.4; 0 0 0; ];

ttb=1000*[0:1/intergrationTime:1/intergrationTime*(size(statesBind,2))]; 
rastergramPlot = subplot(1,1,1,'Parent', figure);
[rasterGramOuput.evi,rasterGramOuput.cev,rasterGramOuput.cia,...
    meanTime2FirstEvent,meanDwellTime,handels.rastergramImg] = ...
    newRastergram(rastergramPlot,statesBind-1,ttb,intergrationTime,...
    align,delete,sortType,sortState,colorArray);


meanDwellTime/1000
meanTime2FirstEvent/1000

idx = statesBind == 1;
TransMatrix2 = modelHMM.fb.transmat;

[states] = genHMMseq(NumTraj, Timesteps, TransMatrix1,ones(NumTraj,1)');
states(idx)=0;
% [tr] = esttrans(states(~idx));
 
rastergramPlot = subplot(1,1,1,'Parent', figure);
[rasterGramOuput.evi,rasterGramOuput.cev,rasterGramOuput.cia,...
    meanTime2FirstEvent,meanDwellTime,handels.rastergramImg] = ...
    newRastergram(rastergramPlot,states,ttb,intergrationTime,...
    align,delete,sortType,sortState,colorArray);

NspotsAn = size(states,1);
NspotsBo=0;

idxFirstEvents = (rasterGramOuput.cia(:,1)==-2);
Time2FirstEvent=rasterGramOuput.cia(idxFirstEvents,5);
idxDwellEvents = (rasterGramOuput.cia(:,1)==1);
DwellTime=rasterGramOuput.cia(idxDwellEvents,5);  

N=16; 
BStrap = 0.1:0.1:0.8;
for j=1:length(BStrap)
    for i=1:N
        idxSelFitst = randperm(size(Time2FirstEvent,1),round(size(Time2FirstEvent,1)*BStrap(j)));
        idxSelDwell = randperm(size(DwellTime,1),round(size(DwellTime,1)*BStrap(j)));
        NspotsBo = NspotsBo +length(idxSelFitst);

        Time2FirstEventTemp=Time2FirstEvent(idxSelFitst);
        time_empty=max(Time2FirstEventTemp)*ones(length(Time2FirstEventTemp),1);
        t=(Time2FirstEventTemp~=time_empty);
        time21st=Time2FirstEventTemp(t);
        meanTime2FirstEventBS(i,j) = mean(time21st);

        meanDwellTimeBS(i,j) = mean(DwellTime(idxSelDwell));

        title(sprintf('spot ID = %d, Arrival time = %0.2g [s] \\pm %0.2g, Dwell time = %0.2g [s] \\pm%0.2g, At %d of %d,\n Analyzed %d Spots, Bootstrapped %d Spots',1,mean(meanTime2FirstEventBS(1:i,j))/1000,...
        std(meanTime2FirstEventBS(1:i,j))/1000,mean(meanDwellTimeBS(1:i,j))/1000,std(meanDwellTimeBS(1:i,j))/1000,i,N,NspotsAn,NspotsBo),'Parent', rastergramPlot)
        meanDwellTimeBS2(i,j) = mean(meanDwellTimeBS(1:i,j));
        meanTime2FirstEventBS2(i,j) = mean(meanTime2FirstEventBS(1:i,j));

        drawnow;
    end
end

%%
figure
m = toff; 
errorbar(mean(meanDwellTimeBS2/1000),std(meanDwellTimeBS2/1000))
hline(m)
hold on
m = ton;
errorbar(mean(meanTime2FirstEventBS2/1000),std(meanTime2FirstEventBS2/1000))
hline(m)

%%

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

prior.m = [200 200 200 200; 200 300 400 500];
sw = load('datatestf2')
mu = sw.mu;
Sigma = sw.Sigma/3;

for i = 1:k
    idx = z==i;
%     Sigma(:,:,i) = iwishrnd(prior.W0(:,:,i),d+1).*eye(2,2); % invpd(wishrnd(W0,v0));
%     mu(:,i) = gaussRnd(prior.m(:,i),prior.kappa(i)*Sigma(:,:,i));
    X(:,idx) = gaussRnd(mu(:,i),Sigma(:,:,i),sum(idx(:)));
end
model.mu = mu;
model.Sigma = Sigma;
model.w = w;
% 
% load('datatestf')
% figure
% for i = 1:k
%     idx = z==i;
%     scatter(X(1,idx),X(2,idx))
%     hold on;
% end
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
% idx = randperm(size(X,2),round(size(X,2)*0.8));
% indexes{i} = idx;
% mu
% 

%%   
%  load('datatest4')
% prior.M  = modelHMM.modelPrior.M(1,1,1);
% prior.alphaA = 1;
% prior.alpha = 1;
% prior.kappa = mean(prior.kappa);
% 
% prior.alpha = 1; 
% 
% prior.v = ones(4,1).*mean(prior.v); 
% prior.W = mean(prior.W0,3)/100; 
% NTimes = 10; 
% NOrder= 4; 
% 
% N=16; 
% BStrap = 0.01:0.1:0.8;
% figure;

% scatter(b(:),c(:))
prior.m = [200.*ones(5,1) 200.*cumsum(ones(5,1));]'; 
prior.m = [200 200 200 200; 200 400+N 600+N 800+N];
X(repmat(permute(states == 0,[3 1 2]),[2 1 1])) = NaN;
order = 1:4;
if size(X,3) > 1 
    dataIdx = ~isnan(X);
    X2 = reshape(X(dataIdx),[size(X,1) sum(sum(sum(dataIdx)))./size(X,1)]);
end
% prior.m = repmat(mean(mean(X2,2),3),[1  max(order)]);
prior.m = [200 200 200 200; 100 300 600 900];
prior.alpha = 1;
prior.kappa = 1
prior.v = 1.*(size(X,1)+1); 
prior.alphaA = 1/2;
prior.M = 1;
[labelest, modelHMMest, evidenceest, modelGMMest] = bayesHMM(X(:,1:50,:), max(order),3,prior,0.05)


%%
dataIdx = ~isnan(X);

[tr] = esttrans(states(squeeze(dataIdx(1,:,:))))
XnoNaN = reshape(X(dataIdx),[size(X,1) sum(sum(sum(dataIdx)))./size(X,1)]);

b = squeeze(XnoNaN(1,:,:));
c = squeeze(XnoNaN(2,:,:));

k=1;
SB=[];
for i=1:size(dataIdx,2)
    startFrames = find(diff(dataIdx(1,:,i)>0) == 1)+1;
    endFrames = find(diff(dataIdx(1,:,i)>0) == -1);
    
    for j=1:size(startFrames,2)
        if endFrames(min(j,end))-startFrames(min(j,end)) >= 1
            if size(endFrames,2) < j
                temp = states(startFrames(j):end,i);
                if length(unique(temp)) > 1
                 SB{k} = temp';
                 k=k+1;
                end

            else
                temp =states(startFrames(j):endFrames(j),i);
                if length(unique(temp)) > 1
                 SB{k} = temp';
                 k=k+1;
                end
            end
        end
    end
end


% x=x-min(x(:))+1;
ns = max([SB{:}]);
% l = size(x,2);
% n = l-1;
traj = size(SB,2);
p = zeros(ns,ns,traj);
for a=1:traj
    l = size(SB{a},2);
    x = SB{a};
    x=x-min(x(:))+1;
    n = l-1;
    for j=1:l
        for t = 1:n
          p(x(1,t), x(1,t + 1),a) = p(x(1,t), x(1,t + 1),a) + 1;
        end
%         for i = 1:ns
%           p(i, :, a) = p(i, :,a) / (eps*10+sum(sum(p(i, :,a))));
%         end
    end
end
[~,x] = max(cellfun(@numel,SB))
pmean = sum(p,3);
mk_stochastic(pmean)

% [tr] = esttrans(SB)
save('datatestf7')