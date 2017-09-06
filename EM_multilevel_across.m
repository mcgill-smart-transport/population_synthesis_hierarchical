function [pi_weight_g, pi_weight_m, log_phi_g, log_phi_m, wj, x_jik, lik] = ...
    EM_multilevel_across(grpid, grpdata, indgid, inddata, G, M, maxIter,con, tol, withTensor, plot_figure,dispiter,initial)

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


if withTensor
    t = num2cell(inddata,1);
    id = sub2ind(indcate,t{:});
    tten = zeros(indcate);
    for i = 1:length(id)
        tten(id(i)) = tten(id(i))+1;
    end
    lraNTF_opts=struct('NumOfComp',M,'CPalg','CP_HALS_opts.mat','lra_rank',M+1,'maxiter',200,'nlssolver','mu','initmode','cp','lra',true);
    res = lraNTF(tten,lraNTF_opts);
    for att = 1:numatt
        w = col_sum(res.U{att});
        res.lambda = res.lambda .* w';
        log_phi_m{att} = log(bsxfun(@rdivide,res.U{att},w));
    end
else
    for att = 1:numatt
        temp = gamrnd(ones(indcate(att),M),con);
        log_phi_m{att} = log(bsxfun(@rdivide, temp, sum(temp,1)));
    end
end
numgat = size(grpdata,2); % number of attributes in group data
log_phi_g = cell(numgat,1);
if withTensor
    t = num2cell(grpdata,1);
    id = sub2ind(grpcate,t{:});
    tten = zeros(grpcate);
    for i = 1:length(id)
        tten(id(i)) = tten(id(i))+1;
    end
    lraNTF_opts=struct('NumOfComp',G,'CPalg','CP_HALS_opts.mat','lra_rank',G+1,'maxiter',200,'nlssolver','mu','initmode','cp','lra',true);
    resg = lraNTF(tten,lraNTF_opts);
    for att = 1:numgat
        w = col_sum(resg.U{att});
        resg.lambda = resg.lambda .* w';
        log_phi_g{att} = log(bsxfun(@rdivide,resg.U{att},w));
    end
    pi_weight_g = resg.lambda/sum(resg.lambda);
    pi_weight_g = pi_weight_g';
else
    for att = 1:numgat
        temp = gamrnd(ones(grpcate(att),G),con);
        log_phi_g{att} = log(bsxfun(@rdivide, temp, sum(temp,1)));
    end
    pi_weight_g = gamrnd(ones(1,G),con);
    pi_weight_g = pi_weight_g/sum(pi_weight_g);
end

nind = size(inddata,1);
ngrp = size(grpdata,1);

pi_weight_m = gamrnd(ones(1,G,M),con);
pi_weight_m = bsxfun(@rdivide,pi_weight_m,sum(pi_weight_m,3));


if ~isempty(initial)
    pi_weight_g = initial.pi_g;
    pi_weight_m = permute(initial.pi_m,[3,1,2]);
    log_phi_g = initial.log_phi_g;
    log_phi_m = initial.log_phi_m;
end
%%

% EM algorithm
iter = 1;
lik = zeros(1,maxIter);
rel = 100;
st = tic;
while iter < maxIter && rel > tol
    % E-step
    % h_jik
    f_ind = zeros(nind,M);
    for att = 1:numatt
        f_ind = f_ind + log_phi_m{att}(inddata(:,att),:);
    end
    f_ind = reshape(f_ind,[nind,1,M]);
    
    x_jik = bsxfun(@times,repmat(pi_weight_m,[nind,1,1]),exp(f_ind));
    % x_jik = pi_weight_m .* exp(f_ind);
    log_g_ji = log(sum(x_jik,3));
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
    
    second = bsxfun(@times,x_jik,wj(indgid,:));
    pi_weight_m = sum(second,1)+eps;
    pi_weight_m = bsxfun(@rdivide, pi_weight_m, sum(pi_weight_m,3));
    
    % M-step
    temp_s = sum(wj,1)+eps;
    pi_weight_g = temp_s/sum(temp_s);
    for att = 1:numgat
        tt = sum_eik(wj, grpdata(:,att), grpcate(att))+eps;
        log_phi_g{att} = log(bsxfun(@rdivide,tt,sum(tt,1)));
    end
    
    second = squeeze(sum(second,2));
    for att = 1:numatt
        tt = sum_eik(second,inddata(:,att),indcate(att))+eps;
        log_phi_m{att} = log(bsxfun(@rdivide,tt,sum(tt,1)));
    end
    % loglikelihood
    temp = bsxfun(@plus,temp_g + tt0, log(pi_weight_g));
    lik(iter) = sum(logsumexp(temp,2));
    if iter > 2
        rel = abs((lik(iter)-lik(iter-1))/lik(iter-1));
    end
    iter = iter + 1;
    
    
    if plot_figure && mod(iter,dispiter) == 0
        t = toc(st);
        fprintf('Iter %i. Average time: %0.3f s. Likelihood %0.3f \n',iter,t/dispiter,lik(iter-1));
        st = tic;
        subplot(2,2,1);
        plot(bsxfun(@times,exp(log_phi_g{1}),pi_weight_g), 'linewidth',2);
        subplot(2,2,2);
        plot(bsxfun(@times,exp(log_phi_g{2}),pi_weight_g),'linewidth',2);
        subplot(2,2,3);
        plot(squeeze(exp(log_phi_m{1})),'linewidth',2);
        subplot(2,2,4);hold on;
        plot(lik(1:iter-1),'-s','linewidth',2);
        xlim([iter*0.2,iter]);
        drawnow;
    end
    
end

disp(iter);
lik = lik(1:iter-1);
