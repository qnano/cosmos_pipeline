clear all

addpath(genpath('../../../helperfunctions'))

load('testmodel','modelGMM','modelHMM')
% load('traces','traces2_allfr', 'ttb')

%%
NumTraj=1000;
Timesteps=5000;
buffer = round(Timesteps*0.1);
d=2;
intergrationTime = 5;

ton = 16;
toff = 8;
TransMatrix1 = [1-(ton*intergrationTime)^-1 (ton*intergrationTime)^-1;  (toff*intergrationTime)^-1 1-(toff*intergrationTime)^-1];
[statesBind] = genHMMseq(NumTraj, Timesteps, TransMatrix1,ones(NumTraj,1)');
[p,psamp] = esttrans(statesBind);

BStrap = 0.1:0.1:0.9;
Nsamps=16
for n=1:Nsamps
    for j=1:length(BStrap)
    idx = randperm(size(psamp,3),round(size(psamp,3)*BStrap(j)));
    length(idx)
    pboot(:,:,j,n) = mean(psamp(:,:,idx),3);
    end
end

figure
pm = mean(pboot,4);
ps = std(pboot,0,4);
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
prior.kappa = ones(k,1); 
prior.m = [200 200 200 200; 350 500 650 800];
for i = 1:k
    idx = z==i;
    Sigma(:,:,i) = iwishrnd(prior.W0(:,:,i),d+1).*eye(2,2); % invpd(wishrnd(W0,v0));
    mu(:,i) = gaussRnd(prior.m(:,i),prior.kappa(i)*Sigma(:,:,i));
    X(:,idx) = gaussRnd(mu(:,i),Sigma(:,:,i),sum(idx(:)));
end

model.mu = mu;
model.Sigma = Sigma;
model.w = w;

plot(squeeze(X(2,10,:))); hold on; plot(0.5.*max(X(2,10,:))*idx(10,:))

%%
prior.M  = modelHMM.modelPrior.M(1,1,1);
prior.alphaA = 1;
prior.alpha = 1;
prior.kappa = mean(prior.kappa);

prior.alpha = 1; 

 
prior.v = ones(4,1).*mean(prior.v); 
prior.W = 100; 
NTimes = 10; 
NOrder= 4; 

N=16; 
BStrap = 0.1:0.1:0.8;

parfor j=1:length(BStrap)
    disp(['j:' num2str(j) 'of ' num2str(length(BStrap))])
    for i=1:N
        idx = randperm(size(X,2),round(size(X,2)*BStrap(j)));
        [labelest{j,i}, modelHMMest{j,i}, evidenceest{j,i}, modelGMMest{j,i}] = bayesHMM(X(:,idx,:), 4,3)
        
       disp(['i:' num2str(i) 'of ' num2str(N)])
    end
end

save('All','-v7.3')
exit