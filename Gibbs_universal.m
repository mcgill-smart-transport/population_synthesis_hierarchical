function [pi_g, pi_m, phi_g, phi_m, alphas, betas,lik] = ...
    Gibbs_universal(data,ind_index,grp_index,indcate,grpcate,NRun,G,M,a_alpha,b_alpha,a_beta,b_beta,a_hyper,dispiter)


% grpid = 1:max(data(:,1));

[~,d] = unique(data(:,1));
grpdata = data(d,grp_index);
indgid = data(:,1);
inddata = data(:,ind_index);

%% initial parameters
numatt = size(inddata,2); % number of attributes in individual data
log_phi_m = cell(numatt,1);
for att = 1:numatt
    temp = gamrnd(ones(indcate(att),M),1);
    log_phi_m{att} = log(bsxfun(@rdivide, temp, sum(temp,1)));
end

numgat = size(grpdata,2); % number of attributes in group data
log_phi_g = cell(numgat,1);
for att = 1:numgat
    temp = gamrnd(ones(grpcate(att),G),1);
    log_phi_g{att} = log(bsxfun(@rdivide, temp, sum(temp,1)));
end

nind = size(inddata,1);
ngrp = size(grpdata,1);

alpha = 1;
beta = 1;

pi_weight_g = gamrnd(ones(1,G),1);
pi_weight_g = pi_weight_g/sum(pi_weight_g);
pi_weight_m = gamrnd(ones(G,M),1);
pi_weight_m = bsxfun(@rdivide,pi_weight_m,sum(pi_weight_m,2));

pi_g = zeros(G, NRun);
pi_m = zeros(G, M, NRun);
phi_g = cell(numgat,1);
alphas = zeros(NRun,1);
betas = zeros(NRun,1);
lik = zeros(NRun,1);
for att = 1:numgat
    phi_g{att} = zeros(grpcate(att),G, NRun);
end
phi_m = cell(numatt,1);
for att = 1:numatt
    phi_m{att} = zeros(indcate(att), M, NRun);
end

%% Gibbs sampling
st = tic;
for iter = 1:NRun
    
    f_ind = zeros(nind,M);
    for att = 1:numatt
        f_ind = f_ind + log_phi_m{att}(inddata(:,att),:);
    end
    
    t = log(pi_weight_m*exp(f_ind'));
    tt0 = zeros(ngrp,G);
    for g = 1:G
        tt0(:,g) = accumarray(indgid, t(g,:));
    end
    
    temp_g = zeros(ngrp,G);
    for att = 1:numgat
        temp_g = temp_g + log_phi_g{att}(grpdata(:,att),:);
    end
    
    temp = bsxfun(@plus,temp_g + tt0, log(pi_weight_g));
    % loglikelihood
    lik(iter) = sum(logsumexp(temp,2));
    
    wj0 = exp(bsxfun(@minus,temp,max(temp,[],2)));
    wj = bsxfun(@rdivide, wj0, sum(wj0,2));
    
    sample_G = mnrnd_new(wj);
    
    temp = sample_G(indgid,:)*pi_weight_m;
    temp = log(temp) + f_ind;
    wi0 = exp(bsxfun(@minus,temp,max(temp,[],2)));
    wi = bsxfun(@rdivide, wi0, sum(wi0,2));
    sample_M = mnrnd_new(wi);
    
    % sample u
    u = zeros(1,G);
    u(G) = 1;
    w = sum(sample_G,1);
    a = 1+w(1:end-1);
    b = alpha+cumsum(w(2:end),'reverse');
    u(1:end-1) = betarnd(a,b);
    t = [1,cumprod(1-u)];
    pi_weight_g = u.*t(1:end-1);
    pi_g(:,iter) = pi_weight_g;
       
    % sample v
    v = zeros(G,M);
    v(:,M) = 1;
    w = sample_G(indgid,:)'*sample_M;
    a = 1+w(:,1:end-1);
    b = beta + cumsum(w(:,2:end),2,'reverse');
    v(:,1:end-1) = betarnd(a,b);
    t = [ones(G,1),cumprod(1-v,2)];
    pi_weight_m = v.*t(:,1:end-1);
    pi_m(:,:,iter) = pi_weight_m;
    
    % sample phi_g
    for att = 1:numgat
        tt = histc(sample_G.*grpdata(:,att),1:grpcate(att));
        tt = gamrnd(tt + a_hyper,1);
        log_phi_g{att} = log(bsxfun(@rdivide,tt,sum(tt,1)));
        phi_g{att}(:,:,iter) = exp(log_phi_g{att});
    end
    
    % sample phi_m
    for att = 1:numatt
        tt = histc(sample_M.*inddata(:,att),1:indcate(att));
        tt = gamrnd(tt + a_hyper,1);
        log_phi_m{att} = log(bsxfun(@rdivide,tt,sum(tt,1)));
        phi_m{att}(:,:,iter) = exp(log_phi_m{att});
    end
    
    % sample alpha
    alpha = gamrnd(a_alpha+G-1, 1/(b_alpha-sum(log(1-u(1:end-1)))));
    alphas(iter) = alpha;
    
    % sample beta
    beta = gamrnd(a_beta+M-1, 1/(b_beta-sum(sum(log(1-v(:,1:end-1))))));
    betas(iter) = beta;
    
    if mod(iter,dispiter)==0
        t = toc(st);
        fprintf('Iter %i. Average time: %0.3f s\n',iter,t/dispiter);
        st = tic;
        subplot(2,2,1);
        [xx,idx] = sort(sum(sample_G,1),'descend');
        plot(xx);
        xlim([0,G]);
        subplot(2,2,2);
        [xx,idy] = sort(sum(sample_M,1),'descend');
        plot(xx);
        xlim([0,M]);
        subplot(2,2,3);
        imagesc(pi_weight_m(idx,idy));
        subplot(2,2,4);
        plot(lik(1:iter));
        xlim([iter*0.2,iter]);
        drawnow;
    end
end