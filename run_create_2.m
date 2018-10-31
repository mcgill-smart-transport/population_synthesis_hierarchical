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

%%
load('res_12_14.mat');
plot(bic);
x = 1; % find the best solution (minimum bic)
pi_m = pi_m{x};
pi_g = pi_g{x};
log_phi_g = log_phi_g{x};
log_phi_m = log_phi_m{x};
wj = wj{x};
xjik = xjik{x};
lik = lik{x};

%% test random generate
close all;
% 10 times of the original PUMS
nhousehold = length(grpid)*10;
num = mnrnd(nhousehold,pi_g);
type = cell(G,1);
for g = 1:G
    type{g} = ones(num(g),1)*g;
end
type = cell2mat(type);

prob = ktensor(pi_g',exp(log_phi_g{1}),exp(log_phi_g{2}),exp(log_phi_g{3}),exp(log_phi_g{4}));
prob_g = cell(G,1);
for g = 1:G
    temp = extract(prob,g);
    temp.lambda = 1;
    prob_g{g} = double(temp);
end

housedata = zeros(nhousehold,length(grplvl));

for g = 1:G
    id = find(type == g);
    thouse = length(id);
    cump = repmat([0;cumsum(prob_g{g}(:))]',thouse,1);
    uf = rand(thouse,1);
    z = zeros(thouse,1);
    for f = 1: numel(prob_g{g})
        z(uf > cump(:,f)) = f;
    end
    ind = cell(1,length(grplvl));
    [ind{:}] = ind2sub(grpcate,z);
    housedata(id,:) = [ind{:}];
end


cate = max(housedata);
temp = num2cell(housedata,1);
ind = sub2ind(cate,temp{:});
counts_syn = zeros(cate);
for i = 1:size(ind,1)
    counts_syn(ind(i)) = counts_syn(ind(i))+1;
end

numAtt = size(housedata,2);
temp = reshape(counts_syn,[],1);
temp = sum(temp,2);
temp = temp/sum(temp);
temp = reshape(temp,cate);
crav = zeros(numAtt,numAtt);
for idxi = 1:numAtt
    for idxj = 1:numAtt
        if idxi ~= idxj
            w = temp;
            for x = 1:numAtt
                if x ~= idxi && x~= idxj, w = sum(w,x); end
            end
            w = squeeze(w);
            w = (w+1e-6)/sum(sum(w));
            w1 = sum(w,1);
            w2 = sum(w,2);
            prd = w2*w1;
            crav(idxi,idxj) = sqrt(sum(sum((w-prd).^2./prd))/min(cate(idxi)-1,cate(idxj)-1));
        end
    end
end
figure; imagesc(crav); crav1 = crav; caxis([0,0.4]);


temp = num2cell(grpdata,1);
ind = sub2ind(cate,temp{:});
counts_raw = zeros(cate);
for i = 1:size(ind,1)
    counts_raw(ind(i)) = counts_raw(ind(i))+1;
end

numAtt = size(housedata,2);
temp = reshape(counts_raw,[],1);
temp = sum(temp,2);
temp = temp/sum(temp);
temp = reshape(temp,cate);
crav = zeros(numAtt,numAtt);
for idxi = 1:numAtt
    for idxj = 1:numAtt
        if idxi ~= idxj
            w = temp;
            for x = 1:numAtt
                if x ~= idxi && x~= idxj, w = sum(w,x); end
            end
            w = squeeze(w);
            w = (w+1e-6)/sum(sum(w));
            w1 = sum(w,1);
            w2 = sum(w,2);
            prd = w2*w1;
            crav(idxi,idxj) = sqrt(sum(sum((w-prd).^2./prd))/min(cate(idxi)-1,cate(idxj)-1));
        end
    end
end
figure; imagesc(crav); crav1 = crav;  caxis([0,0.4])

plot(counts_syn(:),counts_raw(:),'s'); hold on;
plot([0,max(counts_syn(:))],[0,max(counts_syn(:))*sum(counts_raw(:))/sum(counts_syn(:))]);


% generate individuals
clc;
total = sum(housedata(:,4));

persondata = zeros(total,5);
pphouse = housedata(:,4);
% sum(sum(bsxfun(@times,exp(log_phi_m{1}),pi_weight_m),2),3)
prob_m = cell(G,1);
for g = 1:G
    %     prob = ktensor(squeeze(pi_m(1,g,:)),exp(squeeze(log_phi_m{1}(:,g,:))),...
    %         exp(squeeze(log_phi_m{2}(:,g,:))),exp(squeeze(log_phi_m{3}(:,g,:))),...
    %         exp(squeeze(log_phi_m{4}(:,g,:))),exp(squeeze(log_phi_m{5}(:,g,:))));
    prob = ktensor(squeeze(pi_m(1,g,:)),exp(squeeze(log_phi_m{1})),...
        exp(squeeze(log_phi_m{2})),exp(squeeze(log_phi_m{3})),...
        exp(squeeze(log_phi_m{4})),exp(squeeze(log_phi_m{5})));
    prob_m{g} = double(prob);
    sum(prob_m{g}(:))
end


idx = 1;
type2 = cell(nhousehold,1);
for i = 1:nhousehold
    type2{i} = ones(housedata(i,4),1)*type(i);
end
type2 = cell2mat(type2);

x = zeros(14,1);
y = hist(type2,1:G);
for g = 1:G
    x = x+ y(g)*sum(sum(sum(sum(prob_m{g},2),3),4),5);
end
plot(x/sum(x)); hold on;
att = 1;
[y,x] = hist(inddata(:,att),1:max(inddata(:,att))); y = y/sum(y);
plot(x,y,'-ro','linewidth',2); hold on;


for i = 1:G
    id = find(type2 == i);
    len = length(id);
    cump = repmat([0,cumsum(prob_m{i}(:))'],len,1);
    uf = rand(len,1);
    z = zeros(len,1);
    for f = 1: numel(prob_m{i})
        z(uf > cump(:,f)) = f;
    end
    ind = cell(1,length(indcate));
    [ind{:}] = ind2sub(indcate,z);
    persondata(id,:) = [ind{:}];
    %w = [ind{:}];
    %[y,x] = hist(w(:,1),1:14); y = y/sum(y);
    %plot(y); hold on;ma
    %plot(sum(sum(sum(sum(prob_m{i},2),3),4),5),'r');
end


figure;
for att = 1:length(grplvl+1)
    subplot(2,2,att);
    [y,x] = hist(housedata(:,att),1:max(housedata(:,att))); y = y/sum(y);
    plot(x,y,'-s','linewidth',2); hold on;
    [y,x] = hist(grpdata(:,att),1:max(grpdata(:,att))); y = y/sum(y);
    plot(x,y,'-ro','linewidth',2); hold on;
end

figure;
for att = 1:length(indlvl)
    subplot(3,2,att);
    [y,x] = hist(persondata(:,att),1:max(persondata(:,att)));  y = y/sum(y);
    plot(x,y,'-bs','linewidth',2); hold on;
    [y,x] = hist(inddata(:,att),1:max(inddata(:,att))); y = y/sum(y);
    plot(x,y,'-ro','linewidth',2); hold on;
end


type2 = cell(nhousehold,1);
for i = 1:nhousehold
    type2{i} = i*ones(housedata(i,4),1);
end
type2 = cell2mat(type2);

housedata2 = [(1:size(housedata,1))',housedata];
res = [housedata2(type2,1:5),persondata];
res = [res(:,1),(1:size(res,1))',res(:,2:end)];
res(:,6) = [];
%csvwrite('res.csv',res);


grpdata_syn = housedata;
inddata_syn = persondata;
indgid_syn = type2;
%%
cate = max(persondata);
temp = num2cell(persondata,1);
ind = sub2ind(cate,temp{:});
counts_syn = zeros(cate);
for i = 1:size(ind,1)
    counts_syn(ind(i)) = counts_syn(ind(i))+1;
end

numAtt = size(persondata,2);
temp = reshape(counts_syn,[],1);
temp = sum(temp,2);
temp = temp/sum(temp);
temp = reshape(temp,cate);
crav = zeros(numAtt,numAtt);
for idxi = 1:numAtt
    for idxj = 1:numAtt
        if idxi ~= idxj
            x = 1:numAtt;
            w = temp;
            for x = 1:numAtt
                if x ~= idxi && x~= idxj, w = sum(w,x); end
            end
            w = squeeze(w);
            w = (w+1e-6)/sum(sum(w));
            w1 = sum(w,1);
            w2 = sum(w,2);
            prd = w2*w1;
            crav(idxi,idxj) = sqrt(sum(sum((w-prd).^2./prd))/min(cate(idxi)-1,cate(idxj)-1));
        end
    end
end
figure; imagesc(crav); caxis([0,0.4]);



temp = num2cell(inddata,1);
ind = sub2ind(cate,temp{:});
counts_raw = zeros(cate);
for i = 1:size(ind,1)
    counts_raw(ind(i)) = counts_raw(ind(i))+1;
end

numAtt = size(inddata,2);
temp = reshape(counts_raw,[],1);
temp = sum(temp,2);
temp = temp/sum(temp);
temp = reshape(temp,cate);
crav = zeros(numAtt,numAtt);
for idxi = 1:numAtt
    for idxj = 1:numAtt
        if idxi ~= idxj
            x = 1:numAtt;
            w = temp;
            for x = 1:numAtt
                if x ~= idxi && x~= idxj, w = sum(w,x); end
            end
            w = squeeze(w);
            w = (w+1e-6)/sum(sum(w));
            w1 = sum(w,1);
            w2 = sum(w,2);
            prd = w2*w1;
            crav(idxi,idxj) = sqrt(sum(sum((w-prd).^2./prd))/min(cate(idxi)-1,cate(idxj)-1));
        end
    end
end
figure; imagesc(crav); caxis([0,0.4]);

%%
temp = [grpdata(indgid,4),inddata(:,1)];
temp = temp(temp(:,1)==2,:);
temp = abs(temp(2:2:end,2)-temp(1:2:end,2));
hist(temp,0:14);

temp = [housedata(type2,4),persondata(:,1)];
temp = temp(temp(:,1)==2,:);
temp = abs(temp(2:2:end,2)-temp(1:2:end,2));
figure;hist(temp,0:14);

%%
temp = [grpdata(indgid,4),inddata(:,2)];
temp = temp(temp(:,1)==2,:);
temp = [temp(2:2:end,2),temp(1:2:end,2)];
temp = [temp;temp(:,[2,1])];
imagesc(hist3(temp,{1:2,1:2}))

temp = [housedata(type2,4),persondata(:,2)];
temp = temp(temp(:,1)==2,:);
temp = [temp(2:2:end,2),temp(1:2:end,2)];
temp = [temp;temp(:,[2,1])];
figure;imagesc(hist3(temp,{1:2,1:2}));