function [label, model, L] = multiGaussVb(data, m, prior)
% Variational Bayesian inference for HMM with multivariate Gaussian observations.
% Input: 
%   X: d x n data matrix
%   m: k (1 x 1) or label (1 x n, 1<=label(i)<=k) or model structure
% Output:
%   label: 1 x n cluster label
%   model: trained model structure
% Reference: Pattern Recognition and Machine Learning by Christopher M. Bishop (P.474)
% Code is build upon the original VB for mixture of gaussian by Mo Chen (sth4nth@gmail.com).


fprintf('Variational Bayesian Gaussian mixture HMM: running ... \n');
% data = permute(data,[1 3 2]);
% if size(X,3) > 1 
    dataIdx = ~isnan(data);
%     X = reshape(data(dataIdx),[size(data,1) sum(sum(sum(dataIdx)))./size(data,1)]);
% end
X=data;
[num]=size(data,3);

[d,n] = size(dataIdx);
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

model = init(X,dataIdx,m,num,prior);
model = expect(data,dataIdx,model);
model = maximize(X,dataIdx,model,prior);

for iter = 2:maxiter
    model = expect(data,dataIdx,model);
    model = maximize(X,dataIdx,model,prior);
    L(iter) =  bound(data,model,prior); %model.fb.current_loglik/n;
    if abs(L(iter)-L(iter-1)) < tol*abs(L(iter)); break; end
end

L = L(2:iter);
label = zeros(1,sum(sum(sum(dataIdx)))./size(data,1));
[~,label(:)] = max(model.R,[],2);
[~,~,label(:)] = unique(label);
a = reshape(X(dataIdx),[size(X,1) sum(sum(sum(dataIdx)))./size(X,1)]);
figure
scatter(a(1,label==1),a(2,label==1))
hold on
scatter(a(1,label==2),a(2,label==2))

function model = init(X,dataIdx, m,num,prior)
% dataIdx = ~isnan(X);
X = reshape(X(dataIdx),[size(X,1) sum(sum(sum(dataIdx)))./size(X,1)]);
[~, model] = mixGaussVb(X,m,prior);
model.fb.prior = normalise(rand(m,1));
model.fb.transmat = mk_stochastic(ones(m,m));
model.num=num;

% Done
function model = maximize(X,dataIdx, model, prior)
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
m = bsxfun(@plus,kappa0*m0,reshape(X(dataIdx),[size(X,1) sum(sum(sum(dataIdx)))./size(X,1)])*R);
m = bsxfun(@times,m,1./kappa); % 10.61

[d,k] = size(m);
U = zeros(d,d,k); 
logW = zeros(1,k);
r = sqrt(R');
for i = 1:k     
    
    Xm = bsxfun(@times,reshape(X(dataIdx),[size(X,1) sum(sum(sum(dataIdx)))./size(X,1)]),r(i,:));
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
if size(model.fb.prior,1) < size(model.fb.prior,1)
model.fb.prior=    model.fb.prior';
end

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
alpha0 = prior.alpha;
kappa0 = prior.kappa;
v0 = prior.v;
logW0 = prior.logW;
alpha = model.alpha; 
kappa = model.kappa; 
v = model.v;         
logW = model.logW;
R = model.R;
logR = model.logR;
logRho = model.logRho;
[d,n] = size(data);
k = size(R,2);


Elogpi = psi(0,alpha)-psi(0,sum(alpha));
Zz = dot(R(:),logR(:))-dot(R(:),logR(:));

Eppi =  gammaln(k*alpha0)-k*gammaln(alpha0) + (alpha0-1).*sum(Elogpi);
Eqpi =  gammaln(sum(alpha))-sum(gammaln(alpha)) + sum((alpha-1).*Elogpi);
ElogLambda = sum(psi(0,0.5*bsxfun(@minus,v+1,(1:d)')),1)+d*log(2)+logW; % 10.65

Epqmulambda = sum(-0.5*v.*(logW+d*log(2))-logMvGamma(0.5*v,d))-sum(v*d/2)+...
    sum(0.5*(v-d).*ElogLambda)+0.5*sum(log(kappa))-d/2*log(2*pi)-d/2;

Eppmulambda = k.*(-0.5*v0*(logW0+d*log(2))-logMvGamma(0.5*v0,d))-k*v0*d/2+...
    sum(0.5*(v0-d).*ElogLambda)+0.5*k*log(kappa0)-d/2*log(2*pi)-d/2;

L = sum(Zz+Eppi-Eqpi+Epqmulambda-Eppmulambda);
