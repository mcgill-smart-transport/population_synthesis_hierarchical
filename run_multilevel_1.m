% -- Lijun Sun -- %
% -- last modified: September 2, 2017 -- %

clearvars;
clc;
data = csvread('data_raw_person_mat.csv');
% data = csvread('sample.txt');
% data is organized as a merged table of both household attributes and individual attributes
% each row is an individual
% column 1 - household ID; column 2 - individual ID
% column [3,4,5,6] household attributes; column [7,8,9,11,12] individual attributes

data = sortrows(data,[1,2]);
data = bsxfun(@minus,data,min(data)-1);

indgid = data(:,1);
grpid = 1:max(indgid);
numph = int_hist(indgid,max(grpid)); % number of ind per group
data = [data(:,1:5),numph(indgid)',data(:,6:11)];

grplvl = [3,4,5,6];
indlvl = [7,8,9,11,12];


[~,id] = unique(indgid);
grpdata = data(id,grplvl);
inddata = data(:,indlvl);

grpcate = max(grpdata);
indcate = max(inddata);

%% Gibbs
% NRun = 10000;
% G = 20;
% M = 40;
% a_alpha = 0.25; % hyperparameter
% b_alpha = 0.25; % hyperparameter
% a_beta = 0.25; % hyperparameter
% b_beta = 0.25; % hyperparameter
% a_hyper = 1; % hyperparameter
% dispiter = 100;
% 
% [pi_g, pi_m, phi_g, phi_m, alpha, beta,lik] = Gibbs_universal(data,indlvl,grplvl,indcate,grpcate,NRun, G, M,...
%     a_alpha,b_alpha,a_beta,b_beta,a_hyper,dispiter);
% 
% save(strcat('res_gibbs_',num2str(G),'_',num2str(M),'.mat'),'pi_g','pi_m','alpha','beta','phi_g','phi_m','G','lik','M','-v7.3');

%% check performance of the EM code
close all;
profile on;
G = 12;
M = 16;
NRun = 1000;
dispiter = 10;
% [pi_g,pi_m,log_phi_g,log_phi_m,wj,lik] = EM_multilevel_across(grpid', grpdata, indgid, inddata, G, M, ...
%    NRun,10,1e-10,0,1,dispiter,[]);
EM_multilevel_across(grpid', grpdata, indgid, inddata, G, M, ...
    NRun,10,1e-10,0,1,dispiter,[]);
profile viewer;

%% model selection (get the one with minimal BIC value)
% G = 1:25
% M = 6:25
load('bic.mat');
imagesc(mean(bic,3));
imagesc(mean(bic,3));
[a,b,c] = ind2sub(size(bic),find(bic==min(min(min(bic)))));

% optimal G=8, M=16

%% one with the best BIC value
G = 12;
M = 14;
NRUN = 20;
bic = zeros(NRUN,1);
pi_g = cell(NRUN,1);
pi_m = cell(NRUN,1);
log_phi_g = cell(NRUN,1);
log_phi_m = cell(NRUN,1);
wj = cell(NRUN,1);
xjik = cell(NRUN,1);
lik = cell(NRUN,1);
maxIter = 4000;

parfor k = 1:NRUN
    [pi_g{k},pi_m{k},log_phi_g{k},log_phi_m{k},wj{k},xjik{k},lik{k}] = EM_multilevel_across(grpid', grpdata, indgid, inddata, G, M, maxIter,1000000,1e-10,0,1, 100,[]);
    npara = G-1 + G*(M-1) + sum(grpcate-1)*G + sum(indcate-1)*M;
    bic(k) = -2*max(lik{k}) + log(size(inddata,1))*npara;
end
save('res_12_14.mat','bic','pi_g','pi_m','log_phi_g','log_phi_m','wj','xjik','lik','G','M');
%%


