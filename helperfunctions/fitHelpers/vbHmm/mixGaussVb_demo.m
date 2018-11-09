%% Variational Bayesian for Gaussian Mixture Model
close all; clear;
d = 2;
k = 3;
n = 2000;
[X,z] = mixGaussRnd(d,k,n);
plotClass(X,z);

m = floor(n/2);
X1 = X(:,1:m);
X2 = X(:,(m+1):end);
%VB fitting
[y1, model, L] = mixGaussVb(X1,10);
figure;
plotClass(X1,y1);
figure;
plot(L)
all(diff(L)-1e-3<=0)

% Predict testing data
[y2, R] = mixGaussVbPred(model,X2);
figure;
plotClass(X2,y2);

%%

clear all
close all

addpath(genpath('../../helperfunctions/'))


T=1500;
nex = 1; %tracks
M = 1; %g muxture of output
Q = 2; %states


TRANS=[0.8 0.1; 0.1 0.8];
EMIS = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6;...
7/12, 1/12, 1/12, 1/12, 1/12, 1/12];
numTraj=nex;
for i=1:numTraj
    [~,hmmStates(i,:)] = hmmgenerate(T,TRANS,EMIS);
end

% v = ss output
%std output ss
mu1a = 10;
R1a = 50;
z1a = repmat(mu1a,numTraj,T) + randn(numTraj,T)*R1a;

mu1b = 0.1;
R1b = 100;
z1b = repmat(mu1b,numTraj,T) + randn(numTraj,T)*R1b;


mu2a = 600;
R2a = 200;
z2a = repmat(mu2a,numTraj,T) + randn(numTraj,T)*R2a; 


mu2b = 1;
R2b = 100;
z2b = repmat(mu2b,numTraj,T) + randn(numTraj,T)*R2b;

za = (hmmStates == 1).*z1a+(hmmStates == 2).*z2a
zb = (hmmStates == 1).*z1b+(hmmStates == 2).*z2b
seq=[];
seq(:,:,1) = za;
seq(:,:,2) = zb;

data = permute(seq,[3 2 1]);

figure;
% plotClass(data,hmmStates);

% [y1, model, L] = multiGaussVb(data,2);
% [y1, model, L] = mixGaussVbHMM(data, 2)

plotClass(data,y1);
% 
% for ex=1:size(data,3)
%     B = mixgauss_prob(data(:,:,ex), model.m, model.U, [1; 1]);
%     [path(ex,:)] = viterbi_path([0.5;0.5], TRANS, B);
% end
% traces2_allfr_dark=(path-1);
% 
% sum(hmmStates==path)./length(path)

figure;
plot(hmmStates)
hold on
plot(path)

figure
plot(L)




%%
clear all
close all

addpath(genpath('../../helperfunctions/'))

T=100;
nex = 10; %tracks
M = 1; %g muxture of output
Q = 2; %states



TRANS=[0.8 0.1 0.1; 0.25 0.5 0.25; 0.1 0.1 0.8];
EMIS = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6;...
    1/6, 1/6, 1/6, 1/6, 1/6, 1/6;
7/12, 1/12, 1/12, 1/12, 1/12, 1/12];
numTraj=nex;

d = 2;
k = 1;
nc = 2;

[X1,z1] = mixGaussRnd(d,k,T*nex);
[X2,z2] = mixGaussRnd(d,k,T*nex);
[X3,z2] = mixGaussRnd(d,k,T*nex);
% X2=X2+X3+2.5;
X1=reshape(X1,[2 T nex]);
X2=reshape(X2,[2 T nex])+5;
X3=reshape(X3,[2 T nex])+10;
for i=1:numTraj
    [~,hmmStates(i,:,:)] = hmmgenerate(T,TRANS,EMIS);
    data(:,:,i) = (repmat(squeeze(hmmStates(i,:,:))',[d 1])==1).*X1(:,:,i) + (repmat(squeeze(hmmStates(i,:,:))',[d 1])==2).*X2(:,:,i)+ (repmat(squeeze(hmmStates(i,:,:))',[d 1])==3).*X3(:,:,i);
end
%%

[y1, modelTemp3,L1] = mixGaussVbold(reshape(data,[size(data,1) size(data,2)*size(data,3)]), 3);
% [y1, model, L] = mixGaussVb(reshape(data,[size(data,1) size(data,2)*size(data,3)]), 3);

prior.alpha = 1;
prior.alphaA = ones(3);
prior.kappa = 0;
prior.m = [0 5 10;0 5 10];
prior.v = (d+1);
prior.M = 0.01.*repmat(eye(size(data,1)),[1 1 3]);   % M = inv(W)
[y2, model, L] = multiGaussVb(data, 3)
% [v,i] = sort(sum(model.fb.transmat(:,1:nc),2)/nc);


figure
plotClass(data,reshape(permute(hmmStates,[2 3 1]),[1 numel(hmmStates)]));
figure
plotClass(reshape(data,[size(data,1) size(data,2)*size(data,3)]),y1);
figure
plotClass(reshape(data,[size(data,1) size(data,2)*size(data,3)]),y2);


figure
subplot(1,3,1)
scatter(data(1,permute(hmmStates==1,[3 1 2])),data(2,permute(hmmStates==1,[3 1 2])),'or')
subplot(1,3,2)
scatter(data(1,permute(hmmStates==2,[3 1 2])),data(2,permute(hmmStates==2,[3 1 2])),'ob')
subplot(1,3,3)
scatter(data(1,permute(hmmStates==3,[3 1 2])),data(2,permute(hmmStates==3,[3 1 2])),'og')

title('states')

y1=reshape(y1,[1 T nex]);
y2=reshape(y2,[1 T nex]);
figure
subplot(1,3,1)
scatter(data(1,permute(y1==1,[1 2 3])),data(2,permute(y1==1,[1 2 3])),'or')
subplot(1,3,2)
scatter(data(1,permute(y1==2,[1 2 3])),data(2,permute(y1==2,[1 2 3])),'ob')
subplot(1,3,3)
scatter(data(1,permute(hmmStates==3,[3 1 2])),data(2,permute(hmmStates==3,[3 1 2])),'og')

title('mixGaussVbHMM')

figure
subplot(1,3,1)
scatter(data(1,permute(y2==1,[1 2 3])),data(2,permute(y2==1,[1 2 3])),'or')
subplot(1,3,2)
scatter(data(1,permute(y2==2,[1 2 3])),data(2,permute(y2==2,[1 2 3])),'ob')
subplot(1,3,3)
scatter(data(1,permute(hmmStates==3,[3 1 2])),data(2,permute(hmmStates==3,[3 1 2])),'og')

title('multiGauss')

% md=1;
% hmm.trans = model.fb.transmat;
% hmm.prior = model.fb.prior;
% hmmpdf=[];
% mdhmm = hmm2mdhmm(hmm,md);

for ex=1:size(data,3)
    [~, b] = mixGaussVbPred(model, data(:,:,ex));
    [path(ex,:)] = hmmviterbi([],mdhmm,b)
end

%     [B,B2,z] = multiMixGaussVbPred(model, data(:,:,ex));
%     if isfield(mdhmm,'md')
%         b = reshape( repmat(b,mdhmm.md,1), size(b,1), size(b,2)*mdhmm.md);
%     end
%     [path(ex,:)] = viterbi_path(mdhmm.prior,mdhmm.trans, b');

sum(sum(squeeze(hmmStates)==path))./numel(path)
% 
figure;
plot(squeeze(hmmStates(1,1,:)))
hold on
plot(squeeze(path(1,:)),'x')

figure
subplot(1,2,1)
plot(L)
subplot(1,2,2)
plot(L1)

% % now combine states till it becomes worse
% deltall=inf;
% run = 0;
% comb = [];
% while (deltall>0)
% 	% and try if combining helps:
% 	run = run+1;
% 	fprintf('Run %d\n\n',run);
% 	[hmm,pth,comb(run,:),deltall] = hmmcombinestates(hmm,data',pth,stopcrit,reg);
% end
% disp('Experiment done.\n');

% figure
% gplotmatrix(reshape(data,[2 T])',[],reshape(permute(hmmStates,[3 1 2]),[T 1]),['c' 'b' 'm' 'g' 'r'],[],[],false);
% 
% N = 2;                 % initially start with ten segments, each having a 1-component MoG:
% MD = 1;                 % minimum duration of 3 frames
% 
% stopcrit.maxiters = 100;     % stopcriterion: max. nr. of iterations
% stopcrit.minllimpr = 1e-3;  % stopcriterion: min. likelihood improvement
% reg = 1e-3;             % regularization inverse cov. matrices
% 
% %[HMM,PTH,COMB] = hmmtimesegment(data',N,MD,crit,reg);
% 
% disp('Initialize HMM and train\n');
% % split it, and train a hidden node on each of the segments:
% hmm = hmminitwithsegments(data',N);
% % store it in an HMM:
% hmm = hmm2mdhmm(hmm,MD);  % minimum state duration
% disp('Find Viterbi path\n');
% pth = hmmviterbi(data',hmm);
% 
% % now combine states till it becomes worse
% deltall=inf;
% run = 0;
% comb = [];
% while (deltall>0)
% 	% and try if combining helps:
% 	run = run+1;
% 	fprintf('Run %d\n\n',run);
% 	[hmm,pth,comb(run,:),deltall] = hmmcombinestates(hmm,data',pth,stopcrit,reg);
% end
% disp('Experiment done.\n');
%%
close all
evidence =[];
evidence2 =[];
[y2, model, L] = multiGaussVb(data, 3)
% data(:,1501,:)=NaN;

X=reshape(data,[size(data,1) size(data,2)*size(data,3)]);
for i=1:5
    for j=1:5
    prior.alpha = 1;
    prior.alphaA = ones(i).*0.5;
    prior.kappa = 1;sqrt(nex);
    prior.m = 3;
    prior.v = d+1;
    prior.M = repmat(eye(size(data,1)),[1 1 i]);   % M = inv(W)

        [y3, modelTemp3,L4] = mixGaussVb(data, i);
%     [y2, model, L5] = multiGaussVb(data, i,prior)
        evidence2(i,j) = L4(end);
%          evidence(i,j) = L5(end);
%         figure
    subplot(1,2,1)  
            plot(L4)
%     subplot(1,2,2)  
%     plot(L5)
        drawnow;
%         pause(1);
%         [y1, modelTemp3,L3] = mixGaussVb(reshape(data,[size(data,1) size(data,2)*size(data,3)]),i);
%         evidence(i,j) = L3(end);
    end
end
figure
subplot(1,2,2)
plot(mean(evidence,2))
subplot(1,2,1)
plot(mean(evidence2,2))
% 
%%
prior.alpha = 1;
prior.alphaA = 0.5;
prior.kappa = 1;
prior.m = repmat([0 5 10 15 20],[2 1]);
prior.v = d+1;
prior.M = 1;   % M = inv(W)

[label, modelHMM, evidence, modelGMM] = bayesHMM(data, 5,10,prior);
mean(evidence,2)
modelHMM.m
%%

[X1,z1] = mixGaussRnd(d,k,T);
[X2,z2] = mixGaussRnd(d,k,T);

data = (repmat(hmmStates,[d 1])==1).*X1 + (repmat(hmmStates,[d 1])==2).*X2;

[~, modelTemp1] = mixGaussVb(X1,k);
[~, modelTemp2] = mixGaussVb(X2,k);
clear model
model.alphaPhi = permute(cat(1,modelTemp1.alpha,modelTemp2.alpha),[2 1]);
model.kappa = permute(cat(1,modelTemp1.kappa,modelTemp2.kappa),[2 1]);
model.m = permute(cat(3,modelTemp1.m,modelTemp2.m),[1 3 2]);
model.v = permute(cat(1,modelTemp1.v,modelTemp2.v),[2 1]);
model.U = permute(cat(4,modelTemp1.U,modelTemp2.U),[1 2 4 3]);
model.logW = permute(cat(1,modelTemp1.logW,modelTemp2.logW),[2 1]);

model.logRho = permute(cat(3,modelTemp1.logRho,modelTemp2.logRho),[1 2 3]);
model.logR = logsumexp(bsxfun(@minus,model.logRho,logsumexp(logsumexp(model.logRho,2),3)),2); % 10.49
model.R = exp(model.logR);
model.mixMat = ones(k,2);
model.fb.prior = normalise(rand(2,1));
model.fb.transmat = mk_stochastic(rand(2,2));

numex=size(data,3);
exp_num_trans=0;
exp_num_visits1=0;
postmix=0;
for ex=1:numex
    obs = data;
    T = size(obs,2);
     mix = mk_stochastic(model.mixMat);

    [model.Bp,B2] = multiMixGaussVbPred(model,obs);

    [alpha, beta, model.fb.gamma,  model.fb.current_loglik, model.fb.xi_summed,gamma2] =...
        fwdback(model.fb.prior, model.fb.transmat, model.Bp', 'obslik2', permute(B2,[2 3 1]), 'mixmat', mix);

    exp_num_trans = exp_num_trans + model.fb.xi_summed; % sum(xi,3);
    exp_num_visits1 = exp_num_visits1 + model.fb.gamma(:,1);

    postmix = postmix + sum(gamma2,3);
end

for ex=1:size(data,4)
    [B,B2,z] = multiMixGaussVbPred(model, data(:,:,:,ex));
    [path(ex,:)] = viterbi_path(ones(nc,1)'.*0.5, model.fb.transmat, B');
end

traces2_allfr_dark=(path-1);

sum(hmmStates==path)./length(path)

figure;
plot(hmmStates)
hold on
plot(path,'x')