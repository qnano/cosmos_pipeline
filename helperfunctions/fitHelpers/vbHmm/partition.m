function [par,mix] = partition(transmat)

    transmat = mk_stochastic(transmat);
    rowsA = 1:nm*nc;

    for i=0:nc-1
        [colNorm(:,i+1)] = sum(transmat(:,i.*nm+1:(i+1)*nm),2);
        [~,v(:,i+1)] = sort(colNorm(:,i+1),'descend');
    end

    a = [];
    for i=1:length(v(:))
        if ~ismember(v(i),a)
            a=cat(1,a,v(i));
        end
        if length(a) == nc*nm
            break;
        end
    end

    for i=0:nc-1
        [par(i+1,:)] = sum(colNorm(a(i.*nm+1:(i+1)*nm),:),1)/nm;
    end

    for i=0:nc-1
        for j=0:nc-1
            [par(i+1,:)] = sum(colNorm(a(i.*nm+1:(i+1)*nm),:),1)/nm;
            mix{i+1,j+1}=transmat(a(i.*nm+1:(i+1)*nm),j.*nm+1:(j+1)*nm);
        end
    end

%     par
%     mix

end