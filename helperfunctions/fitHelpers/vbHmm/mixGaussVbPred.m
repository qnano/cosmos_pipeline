function [z, R] = mixGaussVbPred(model, X)
% Predict label and responsibility for Gaussian mixture model trained by VB.
% Input:
%   X: d x n data matrix
%   model: trained model structure outputed by the EM algirthm
% Output:
%   label: 1 x n cluster label
%   R: k x n responsibility
% Written by Mo Chen (sth4nth@gmail.com).
alpha = model.alpha; % Dirichlet
kappa = model.kappa;   % Gaussian
m = model.m;         % Gasusian
v = model.v;         % Whishart
U = model.U;         % Whishart 
logW = model.logW;
n = size(X,2);
[d] = size(m,1);
k=size(U,3);
EQ = zeros(n,k);
for i = 1:k
    for j=1:size(X,3)
        Q = (U(:,:,i)'\bsxfun(@minus,squeeze(X(:,:,j)),m(:,i)));
        EQ(:,i,j) = d/kappa(i)+v(i)*dot(Q,Q,1);    % 10.64
    end
end
ElogLambda = sum(psi(0,0.5*bsxfun(@minus,v+1,(1:size(v,1))')),1)+d*log(2)+logW; % 10.65
Elogpi = psi(0,alpha)-psi(0,sum(alpha)); % 10.66
for j=1:size(X,3)
    logRho(:,:,j) = -0.5*bsxfun(@minus,squeeze(EQ(:,:,j)),ElogLambda(:,:,j)-d*log(2*pi)); % 10.46
    %logRho(:,:,j) = bsxfun(@plus,logRho(:,:,j),Elogpi(:,:,j));   % 10.46
    logR(:,:,j) = bsxfun(@minus,logRho(:,:,j),logsumexp(logRho(:,:,j),2)); % 10.49
    R(:,:,j) = exp(logR(:,:,j));
%     z = zeros(size(X,3),n);
    [~,z(j,:)] = max(R(:,:,j),[],2);
    [~,~,z(j,:)] = unique(z(j,:));
end

