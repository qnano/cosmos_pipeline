function [sortedModel ] = sortModel(modelHMM,I)
    sortedModel.R = modelHMM.R(:,I);
    sortedModel.alpha = modelHMM.alpha(I);
    sortedModel.kappa = modelHMM.kappa(I);
    sortedModel.m = modelHMM.m(:,I);
    sortedModel.v = modelHMM.v(I);
    sortedModel.U = modelHMM.U(:,:,I);
    sortedModel.logW = modelHMM.logW(I);
    sortedModel.logR = modelHMM.logR(:,I);
    sortedModel.logRho = modelHMM.logRho(:,I);
    sortedModel.num = modelHMM.num;
    sortedModel.alphaA = modelHMM.alphaA(I,I);

    sortedModel.fb.prior = modelHMM.fb.prior(I);
    sortedModel.fb.current_loglik = modelHMM.fb.current_loglik;
    sortedModel.fb.gamma = modelHMM.fb.gamma(I,:);
    sortedModel.fb.obs = modelHMM.fb.obs(I,:);
    sortedModel.fb.exp_num_visits1 = modelHMM.fb.exp_num_visits1(I);
    sortedModel.fb.transmat = modelHMM.fb.transmat(I,I);
    sortedModel.fb.xi_summed = modelHMM.fb.xi_summed(I,I);
    sortedModel.fb.xi = modelHMM.fb.xi(I,I,:);
end