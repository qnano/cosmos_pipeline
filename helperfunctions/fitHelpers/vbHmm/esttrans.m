function [tr] = esttrans(statesBind)
%     x=x-min(x(:))+1;
%     ns = max(x(:));
%     l = size(x,2);
%     n = l-1;
%     traj = size(x,1);
%     p = zeros(ns,ns,traj);
%     for a=1:traj
%         for j=1:l
%             for t = 1:n
%               p(x(a,t), x(a,t + 1),a) = p(x(a,t), x(a,t + 1),a) + 1;
%             end
%             for i = 1:ns
%               p(i, :, a) = p(i, :,a) / (eps*10+sum(sum(p(i, :,a))));
%             end
%         end
%     end
%     pmean = mean(p,3);
k=1;
if ismatrix(statesBind)
    for i=1:size(statesBind,1)
        idxa = (sum(statesBind(i,:),1)>0)';
        starts = find(idxa.*diff([0;idxa]));
        ends = find(abs(idxa.*diff([idxa;0])));
        for j=1:size(starts,1)
            idxa = starts(j):ends(j);
            SB(k) = {squeeze(statesBind(i,idxa))};
            k=k+1;
        end
    end
else
    SB = stateBind;
end
transitions = cellfun(@(x)([x(1:length(x)-1); x(2:length(x))]), SB, 'UniformOutput', false);

alltransitions = cell2mat(transitions)';


[uniqueTransitions, ~, q]=unique(alltransitions,'rows','stable');
v=arrayfun(@(x) sum(q==x),1:size(uniqueTransitions,1))';
p = v/sum(v);
nmax = max(max(uniqueTransitions));
transitionMatrix = sparse(uniqueTransitions(:,1), uniqueTransitions(:,2), p, nmax,nmax);
tr = mk_stochastic(full(transitionMatrix));
