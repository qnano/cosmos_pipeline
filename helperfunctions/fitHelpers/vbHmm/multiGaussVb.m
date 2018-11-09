function [labelposterior, model, labellikelihood,L] = multiGaussVb(data, m, prior)
% Variational Bayesian inference for Gaussian mixture.
% Input: 
%   X: d x n data matrix
%   m: k (1 x 1) or label (1 x n, 1<=label(i)<=k) or model structure
% Output:
%   label: 1 x n cluster label
%   model: trained model structure


fprintf('Variational Bayesian Gaussian mixture HMM: running ... \n');
[num]=size(data,3);
dataIdx = ~isnan(data);
X = reshape(data(dataIdx),[size(data,1) sum(sum(sum(dataIdx)))./size(data,1)]);

[d,n] = size(X);
if nargin < 3
    prior.alpha = 1;
    prior.alphaA = ones(m);
    prior.kappa = 1;
    prior.m = repmat(mean(X,2),[1 m]);
    prior.v = d+1;
    prior.M = repmat(eye(d),[1 1 m]);   % M = inv(W)
end
for i=1:size(prior.M,3)
    prior.logW(i) = -2*sum(log(diag(chol(prior.M(:,:,i)))));
end
prior.logW=mean(prior.logW);

tol = 1e-8;
maxiter = 200;
L = -inf(1,maxiter);

model = init(X,m,num,prior);
model = expect(data,dataIdx,model);
model = maximize(X,model,prior);

for iter = 2:maxiter
    model = expect(data,dataIdx,model);
    model = maximize(X,model,prior);
    L(iter) = model.fb.current_loglik/n;
    if abs(L(iter)-L(iter-1)) < tol*abs(L(iter)); break; end
end

L = L(2:iter);
labelposterior = zeros(1,n);
labellikelihood = zeros(1,n);
[v,I] = sort(mean(model.m,1));
model = sortModel(model,I);

ElogA = psi(0,model.alphaA)-psi(0,sum(model.alphaA,2));
Zz= exp(model.fb.obs.*model.R')+squeeze(sum(model.fb.xi.*repmat(ElogA,[1 1 size(model.fb.xi,3) size(model.fb.xi,4)]),2));
Zzz = Zz./sum(Zz,1);
[~,labelposterior(:)] = max(Zzz,[],1);
[~,~,labelposterior(:)] = unique(labelposterior);

[~,labellikelihood(:)] = max(model.R,[],2);
[~,~,labellikelihood(:)] = unique(labellikelihood);

% a = reshape(data(dataIdx),[size(data,1) sum(sum(sum(dataIdx)))./size(data,1)]);
% figure
% scatter(a(1,label==1),a(2,label==1))
% hold on
% scatter(a(1,label==2),a(2,label==2))

function model = init(X, m,num,prior)
[~, model] = mixGaussVb(X,m,prior);
model.fb.prior = normalise(rand(m,1));
model.fb.transmat = mk_stochastic(rand(m,m));
model.num=num;

% Done
function model = maximize(X, model, prior)
alpha0 = prior.alpha;
alphaA0 = prior.alphaA;
kappa0 = prior.kappa;
m0 = prior.m;
v0 = prior.v;
M0 = prior.M;

R = model.fb.gamma';
RR = model.fb.xi_summed;

nk = sum(R,1); % 10.51
alpha = alpha0+model.fb.exp_num_visits1'; % 10.58
alphaA=alphaA0+RR; %check
kappa = kappa0+nk; % 10.60
v = v0+nk; % 10.63
m = bsxfun(@plus,kappa0*m0,X*R);
m = bsxfun(@times,m,1./kappa); % 10.61

[d,k] = size(m);
U = zeros(d,d,k); 
logW = zeros(1,k);
r = sqrt(R');
for i = 1:k     
    Xm = bsxfun(@times,X,r(i,:));
    M = M0(:,:,i)+Xm*Xm'+kappa0*(m0(:,i)*m0(:,i)')-kappa(i)*(m(:,i)*m(:,i)');     % equivalent to 10.62
    U(:,:,i) = chol(M);
    logW(i) = -2*sum(log(diag(U(:,:,i))));      
end

model.alpha = alpha;
model.alphaA = alphaA;
model.kappa = kappa;
model.m = m;
model.v = v;
model.U = U;
model.logW = logW;
model.fb.transmat = mk_stochastic(exp(psi(0,model.alphaA)-psi(0,sum(model.alphaA,2))));
model.fb.prior= normalise(exp(psi(0,model.alpha)-psi(0,sum(model.alpha))));


function model = expect(data,dataIdx, model)
exp_num_trans=0;
exp_num_visits1=0;
current_loglik=0;

N = sum(sum(sum(dataIdx)));
gamma = nan(size(model.R,2),size(data,2),model.num);
xi = nan(size(model.R,2),size(model.R,2),size(data,2)-1,model.num);
obs = nan(size(model.R,2),size(data,2),model.num);
ntraces=0;
for ex=1:model.num
    idxa = (sum(dataIdx(:,:,ex),1)>0)';
    starts = find(idxa.*diff([0;idxa]));
    ends = find(abs(idxa.*diff([idxa;0])));
    for j=1:size(starts,1)
        idxa = starts(j):ends(j);
        [~, B] = mixGaussVbPred(model, data(:,idxa,ex));
        obs(:,idxa,ex)=permute(B, [2 1]);
        [~, ~, gamma(:,idxa,ex),  current_logliknew, xi_summed,xi(:,:,idxa(1:end-1),ex)]  = fwdback(model.fb.prior, model.fb.transmat, B');
        exp_num_trans = exp_num_trans + xi_summed; 
        exp_num_visits1 = exp_num_visits1 + gamma(:,starts(j),ex);
        current_loglik=current_loglik+current_logliknew;
        ntraces = ntraces+1;
    end
end

model.fb.current_loglik=current_loglik;
model.fb.xi_summed = exp_num_trans;
idxGamma = ~isnan(gamma);
model.fb.gamma=reshape(gamma(idxGamma),[size(gamma,1) sum(sum(sum(idxGamma)))./size(gamma,1)]);
model.fb.obs=reshape(obs(idxGamma),[size(gamma,1) sum(sum(sum(idxGamma)))./size(gamma,1)]);
model.fb.exp_num_visits1=exp_num_visits1/ntraces;
xi(:,:,end+1,:) = 0;
xi = circshift(xi,[0 0 1 ]);
model.fb.xi = reshape(xi(:,:,idxGamma(1,1:end,:)),[size(xi,1) size(xi,2)  sum(sum(idxGamma(1,1:end,:)))]);


function L = bound(data,model,prior)
    
    L=model.fb.current_loglik;
