clear all
close all

BStrap = 0.01:0.1:0.8;
N=64; 
for j=1:length(BStrap)
    ws = load(['HMMloc2' num2str(j)]);
    disp(['j:' num2str(j) 'of ' num2str(length(BStrap))])
    for i=1:N
%        disp(['i:' num2str(i) 'of ' num2str(N)])        
       A{i,j} = ws.modelHMMest{i}.fb.transmat;
       NOrd(i,j) = size(A{i,j},1);
       samples{i,j} = ws.X(:,ws.indexes{i},:);
       idexes{i,j}=ws.indexes{i};
       label{i,j} = ws.labelest{i};
       occupationcy{i,j} = round(sum(ws.modelHMMest{i}.R,1)./size(ws.modelHMMest{i}.R,1)*100)/100;
%        ws2.states(idexes{i,j}',:)      
       
       a = samples{i,j};
       b = a(~isnan(a));
       c = reshape(b,[2 size(b,1)/2]);
       for k=1:max(label{i,j}(:))
           d = c(:,label{i,j}== k);
           SigmaEst(:,:,i,j,k) = cov(d');
           muEst(:,i,j,k) = mean(d');
       end
    end
end

% ws2 = load('nicerdata.mat');
for i=1:max(ws.states(:))
    a = ws.X(:,ws.states == i);
    sigmaTrue(:,:,i) = cov(a');
    muTrue(:,i) = mean(a');
end
sigmaTrue = ws.Sigma;
muTrue = ws.mu;
    
for j=1:length(BStrap)
    ws = load(['HMMloc2qq1' num2str(j)]);
    disp(['j:' num2str(j) 'of ' num2str(length(BStrap))])
    for i=1:N
%        disp(['i:' num2str(i) 'of ' num2str(N)])        
       A{i+N,j} = ws.modelHMMest{i}.fb.transmat;
       NOrd(i+N,j) = size(A{i,j},1);
       samples{i+N,j} = ws.X(:,ws.indexes{i},:);
       label{i+N,j} = ws.labelest{i};
       occupationcy{i+N,j} = mat2str(round(sum(ws.modelHMMest{i}.R,1)./size(ws.modelHMMest{i}.R,1)*100)/100);
       
       a = samples{i+N,j};
       b = a(~isnan(a));
       c = reshape(b,[2 size(b,1)/2]);
       for k=1:max(label{i+N,j}(:))
           d = c(:,label{i+N,j}== k);
           SigmaEst(:,:,i+N,j,k) = cov(d');
           muEst(:,i+N,j,k) = mean(d');
       end
    end
end
N=128
M=64
for j=1:length(BStrap)
    ws = load(['HMMloc4' num2str(j)]);
    disp(['j:' num2str(j) 'of ' num2str(length(BStrap))])
    for i=1:N
%        disp(['i:' num2str(i) 'of ' num2str(N)])        
       A{i+M,j} = ws.modelHMMest{i}.fb.transmat;
       NOrd(i+M,j) = size(A{i,j},1);
       samples{i+M,j} = ws.X(:,ws.indexes{i},:);
       label{i+M,j} = ws.labelest{i};
       occupationcy{i+M,j} = mat2str(round(sum(ws.modelHMMest{i}.R,1)./size(ws.modelHMMest{i}.R,1)*100)/100);
       
       a = samples{i+M,j};
       b = a(~isnan(a));
       c = reshape(b,[2 size(b,1)/2]);
       for k=1:max(label{i+M,j}(:))
           d = c(:,label{i+M,j}== k);
           SigmaEst(:,:,i+M,j,k) = cov(d');
           muEst(:,i+M,j,k) = mean(d');
       end
    end
end
% 
N=128
M=128+64
for j=1:length(BStrap)
    ws = load(['HMMloc3' num2str(j)]);
    disp(['j:' num2str(j) 'of ' num2str(length(BStrap))])
    for i=1:N
%        disp(['i:' num2str(i) 'of ' num2str(N)])        
       A{i+M,j} = ws.modelHMMest{i}.fb.transmat;
       NOrd(i+M,j) = size(A{i,j},1);
       samples{i+M,j} = ws.X(:,ws.indexes{i},:);
       label{i+M,j} = ws.labelest{i};
       occupationcy{i+M,j} = mat2str(round(sum(ws.modelHMMest{i}.R,1)./size(ws.modelHMMest{i}.R,1)*100)/100);
       
       a = samples{i+M,j};
       b = a(~isnan(a));
       c = reshape(b,[2 size(b,1)/2]);
       for k=1:max(label{i+M,j}(:))
           d = c(:,label{i+M,j}== k);
           SigmaEst(:,:,i+M,j,k) = cov(d');
           muEst(:,i+M,j,k) = mean(d');
       end
    end
end

%  mat2str(round(sum(ws.modelHMMest{i}.R,1)./size(ws.modelHMMest{i}.R,1)*100)/100)

for i=1:max(NOrd(:))
    v(i,:) = sum(NOrd' == i,2)/N;
end

%%
a = [220 44 83; 64 78 200; 234 193 30; 64 200 78 ]/255;

H = bar(v', 'stacked', 'BarWidth', 1)
colormap(cool)
% colorSet = [];
% for i = 1:3
%     myColors = a;
%     colorSet = [colorSet myColors];
% %     H(i).FaceColor = 'flat';
%     H(i).FaceColor = myColors(i,:);
% end
% axis tight
% 
% legend

%%


figure
pm = mean(ws.pboot,4);
ps = std(ws.pboot,0,4);

subplot(2,2,1) 
errorbar(squeeze(pm(1,1,:)),squeeze(ps(1,1,:)))
hline(ws.TransMatrix1(1,1))
ylim([0.98 1])
ylabel('A_{1,1} [s^-1]')

subplot(2,2,2)
errorbar(squeeze(pm(1,2,:)),squeeze(ps(1,2,:)))
hline(ws.TransMatrix1(1,2))
ylabel('A_{1,2} [s^-1]')

subplot(2,2,3)
errorbar(squeeze(pm(2,1,:)),squeeze(ps(2,1,:)))
hline(ws.TransMatrix1(2,1))
ylabel('A_{2,1} [s^-1]')

subplot(2,2,4)
errorbar(squeeze(pm(2,2,:)),squeeze(ps(2,2,:)))
hline(ws.TransMatrix1(2,2))
ylabel('A_{2,2} [s^-1]')


% Kullback Leibler divergence (KL) KL(p,q) between two gaussian
% distributions:
clear dKL
for j=1:size(SigmaEst,3)
    for k=1:size(SigmaEst,4)
        sigma1 = squeeze(SigmaEst(:,:,j,k,1:2));
        mu1 = squeeze(muEst(:,j,k,1:2));
        dKL(j,k,:) = sqrt(calculate_distance(sigma1,mu1,sigmaTrue,muTrue));
    end
end

figure
a = ~isnan(sum(dKL(:,:,:),3));

for i=1:size(a,2)
    ss = dKL(a(:,i),i,:);
    MKL(i,:) = squeeze(mean(ss,1));
    SKL(i,:) = squeeze(std(ss,1,1));
end

errorbar(MKL,SKL)


%%
figure;
i=3
subplot(1,3,1)
scatter(ws.X(1,ws.states==1),ws.X(2,ws.states==1))
hold on
scatter(ws.X(1,ws.states==2),ws.X(2,ws.states==2))
subplot(1,3,2)
a = ws.X(:,ws.indexes{i},:);
[y1, model, L] = mixGaussVb(a,2);
dataIdx = ~isnan(a);
b = reshape(a(dataIdx),[size(a,1) sum(sum(sum(dataIdx)))./size(a,1)]);
% c = ws.labelest{i};
c=y1;
scatter(b(1,c==1),b(2,c==1))
hold on
scatter(b(1,c==2),b(2,c==2))

subplot(1,3,3)
[y2, modelHMM,L] = multiGaussVb(a, 2);
dataIdx = ~isnan(a);
b = reshape(a(dataIdx),[size(a,1) sum(sum(sum(dataIdx)))./size(a,1)]);
c = y2;
scatter(b(1,c==1),b(2,c==1))
hold on
scatter(b(1,c==2),b(2,c==2))

figure;

plot(L)

cov(b(:,c==1)')
cov(b(:,c==2)')