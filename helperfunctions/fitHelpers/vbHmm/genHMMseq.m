% genHMMseq generates a Markov chain.
% [states]= genHMMseq(NumTraj, Timesteps, TransMatrix)
% 
% input: 
% NumTraj: the number of different trajectories simulated (integer!)
% Timesteps: the length of the trajectories
% TransMatrix: transition matrix [p11 p12 p13...; p21 p22 ...; p31 p32...; ... 
%    p41 p42..]. Must be a square matrix

% output: 
% states: state sequences. Is a matrix with the dimensions [NumTraj Timesteps]

% Carlas Smith, October 2014


function [states] = genHMMseq(NumTraj, Timesteps, TransMatrix, s0)
% tr must be square

numStates = size(TransMatrix,1);
checkTr = size(TransMatrix,2);
if checkTr ~= numStates;
    error(message('stats:hmmgenerate:BadTransitions'));
end

% calculate cumulative probabilities
trc = cumsum(TransMatrix,2);

% normalize these just in case they don't sum to 1.
trc = trc./repmat(trc(:,end),1,numStates);



% create random sequence for state changes
statechange = rand(NumTraj,Timesteps-1);
states=zeros(NumTraj,Timesteps);

if nargin < 4
    [V,D,W]=eig(TransMatrix );
    [a,b]=find(single(D)==1);    perc= cumsum(W(:,b)./sum(W(:,b)));
    U=rand(NumTraj,1);
    
for i=1:NumTraj
    states(i,1)=find(U(i)<=perc,1);
end
    
else
    states(:,1)=s0;
end

for j=1:NumTraj
    for i=1:Timesteps-1
    states(j,i+1)= find(statechange(j,i)<=trc(states(j,i),:),1);
    end
end

