function [pi_weight_g, pi_weight_m, log_phi_g, log_phi_m, wj, x_jik, lik] = ...
    EM_multilevel(grpid, grpdata, indgid, inddata, G, M, maxIter,con, tol)

% con: alpha
% G: grp latent classes
% M: ind latent classes

% grpid group ID
% indgid group ID in ind data
% total number of missing variables is very large

grpcate = max(grpdata);
indcate = max(inddata);

numatt = size(inddata,2); % number of attributes in individual data
log_phi_m = cell(numatt,1);
for att = 1:numatt
    temp = gamrnd(ones(indcate(att),G,M),con);
    log_phi_m{att} = log(bsxfun(@rdivide, temp, sum(temp,1)));
end

numgat = size(grpdata,2); % number of attributes in group data
log_phi_g = cell(numgat,1);
for att = 1:numgat
    temp = gamrnd(ones(grpcate(att),G),con);
    log_phi_g{att} = log(bsxfun(@rdivide, temp, sum(temp,1)));
end

nind = size(inddata,1);
ngrp = size(grpdata,1);


pi_weight_g = gamrnd(ones(1,G),con);
pi_weight_g = pi_weight_g/sum(pi_weight_g);
pi_weight_m = gamrnd(ones(1,G,M),con);
pi_weight_m = bsxfun(@rdivide,pi_weight_m,sum(pi_weight_m,3));
%%

% EM algorithm
iter = 1;
lik = zeros(1,maxIter);
rel = 100;
while iter < maxIter && rel > tol
    % E-step
    % h_jik
    f_ind = zeros(nind,G,M);
    for att = 1:numatt
        f_ind = f_ind + log_phi_m{att}(inddata(:,att),:,:);
    end
    
    % h_jik
    log_h_jik = f_ind + repmat(log(pi_weight_m),[nind,1,1]);
    
    log_g_ji = logsumexp(log_h_jik,3);
    
    log_x_jik = bsxfun(@minus, log_h_jik, max(log_h_jik,[],3));
    x_jik = exp(log_x_jik);
    x_jik = bsxfun(@rdivide,x_jik,sum(x_jik,3));
    
    tt0 = zeros(ngrp,G);
    for g = 1:G
        tt0(:,g) = accumarray(indgid, log_g_ji(:,g));
    end
    
    temp_g = zeros(ngrp,G);
    for att = 1:numgat
        temp_g = temp_g + log_phi_g{att}(grpdata(:,att),:);
    end
    
    temp = bsxfun(@plus,temp_g + tt0, log(pi_weight_g));
    wj0 = exp(bsxfun(@minus,temp,max(temp,[],2)));
    wj = bsxfun(@rdivide, wj0, sum(wj0,2));
    
    second = times(x_jik,repmat(wj(indgid,:),[1,1,M]));
    pi_weight_m = sum(second,1);
    pi_weight_m = bsxfun(@rdivide, pi_weight_m, sum(pi_weight_m,3));
    
    % M-step
    temp_s = sum(wj,1);
    pi_weight_g = temp_s/sum(temp_s);
    for att = 1:numgat
        tt = sum_eik(wj, grpdata(:,att), grpcate(att));
        log_phi_g{att} = log(bsxfun(@rdivide,tt,sum(tt,1)));
    end
    
    for att = 1:numatt
        for k = 1:M
            tt = sum_eik(second(:,:,k),inddata(:,att),indcate(att));
            log_phi_m{att}(:,:,k) = log(bsxfun(@rdivide,tt,sum(tt,1)));
        end
        %log_phi_m{att}(log_phi_m{att}<=-500) = -100000;
    end
    % loglikelihood
    temp = bsxfun(@plus,temp_g + tt0, log(pi_weight_g));
    lik(iter) = sum(logsumexp(temp,2));
    if iter > 2
        rel = abs((lik(iter)-lik(iter-1))/lik(iter-1));
    end
    iter = iter + 1;
    
    if mod(iter,10) == 0
        subplot(2,2,1);
        plot(bsxfun(@times,exp(log_phi_g{1}),pi_weight_g), 'linewidth',2);
        subplot(2,2,2);
        plot(bsxfun(@times,exp(log_phi_g{2}),pi_weight_g),'linewidth',2);
        subplot(2,2,3);
        %plot(squeeze(exp(log_phi_m{1})),'linewidth',2);
        subplot(2,2,4);
        plot(lik(1:iter-1),'-s','linewidth',2);
        xlim([iter*0.2,iter]);
        drawnow;
        disp([iter,lik(iter-1)]);
    end
end

disp(iter);
lik = lik(1:iter-1);


