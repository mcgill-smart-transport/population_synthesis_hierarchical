% -- Lijun Sun -- %
% -- last modified: September 2, 2017 -- %

h = fig('units','centimeters','width',18,'height',5,'font','helvetica','fontsize',7,'border','off');
tempr = [grpdata(indgid,4),inddata];
tempr = tempr(tempr(:,1)==2,:);
ass1 = [tempr(1:2:end,2:end),tempr(2:2:end,2:end);...
    tempr(2:2:end,2:end),tempr(1:2:end,2:end)];

numAtt = size(ass1,2);
cate = max(ass1);
crav = zeros(numAtt,numAtt);
for idxi = 1:numAtt-1
    for idxj = idxi+1:numAtt
        w = hist3(ass1(:,[idxi,idxj]),{1:cate(idxi),1:cate(idxj)});
        w = (w+1e-6)/sum(w(:));
        w1 = sum(w,1);
        w2 = sum(w,2);
        prd = w2*w1;
        crav(idxi,idxj) = sqrt(sum(sum((w-prd).^2./prd))/min(cate(idxi)-1,cate(idxj)-1));
    end
end
crav1 = crav+crav';
subplot(1,3,1); imagesc(crav1);  caxis([0,0.7]); hold on;
plot([5.5,5.5],[0,11],'-','color','w','linewidth',0.5); hold on;
plot([0,11],[5.5,5.5],'-','color','w','linewidth',0.5); hold on;
set(gca,'xtick',1:10,'xticklabel',{'A','S','I','L','P','A','S','I','L','P'});
set(gca,'ytick',1:10,'yticklabel',{'A','S','I','L','P','A','S','I','L','P'});
xlabel('\bf{Original}','interpreter','latex');

temp = [grpdata_syn(indgid_syn,4),inddata_syn];
temp = temp(temp(:,1)==2,:);
ass2 = [temp(1:2:end,2:end),temp(2:2:end,2:end);...
    temp(2:2:end,2:end),temp(1:2:end,2:end)];
crav = zeros(numAtt,numAtt);
for idxi = 1:numAtt-1
    for idxj = idxi+1:numAtt
        w = hist3(ass2(:,[idxi,idxj]),{1:cate(idxi),1:cate(idxj)});
        w = (w+1e-6)/sum(w(:));
        w1 = sum(w,1);
        w2 = sum(w,2);
        prd = w2*w1;
        crav(idxi,idxj) = sqrt(sum(sum((w-prd).^2./prd))/min(cate(idxi)-1,cate(idxj)-1));
    end
end
crav2 = crav+crav';
subplot(1,3,2); imagesc(crav2);  caxis([0,0.7]); hold on;
plot([5.5,5.5],[0,11],'-','color','w','linewidth',0.5); hold on;
plot([0,11],[5.5,5.5],'-','color','w','linewidth',0.5); hold on;
set(gca,'xtick',1:10,'xticklabel',{'A','S','I','L','P','A','S','I','L','P'});
set(gca,'ytick',1:10,'yticklabel',{'A','S','I','L','P','A','S','I','L','P'});
xlabel('\bf{Synthetic}','interpreter','latex');

subplot(1,3,3); imagesc(abs(crav1-crav2)); caxis([0,0.7]); hold on;
plot([5.5,5.5],[0,11],'-','color','w','linewidth',0.5); hold on;
plot([0,11],[5.5,5.5],'-','color','w','linewidth',0.5); hold on;
set(gca,'xtick',1:10,'xticklabel',{'A','S','I','L','P','A','S','I','L','P'});
set(gca,'ytick',1:10,'yticklabel',{'A','S','I','L','P','A','S','I','L','P'});
xlabel('\bf{Difference}','interpreter','latex');

hh = colorbar('fontsize',6);
h_bar = findobj(gcf,'Tag','Colorbar');
initpos = get(h_bar,'Position');
pos2 = [0.93,0.25,0.02,0.65];
set(h_bar, 'Position',pos2);
set(h_bar,'ticklength',0.015,'LineWidth',0.1);
set(h_bar,'ylim',[0,0.7]);

name = 'fig_compare_two_people.pdf';
export_fig(name,'-pdf','-transparent','-painters','-rgb');
%close all;
system(['start ',name]);


%%
clc;
tempr = [grpdata(indgid,4),inddata];
tempr = tempr(tempr(:,1)==2,:);
% age1, age2, sex1, sex2
temp2r = [tempr(2:2:end,2),tempr(1:2:end,2),tempr(2:2:end,3),tempr(1:2:end,3)];
temp2r = [temp2r;temp2r(:,[2,1,4,3])];
% [max(age), diff(age), sex1+sex2]
temp2 = [max(temp2r(:,[1,2]),[],2), abs(temp2r(:,1)-temp2r(:,2))+1,(temp2r(:,3)-1)*2+temp2r(:,4)];
temp2 = [max(temp2r(:,[1,2]),[],2), abs(temp2r(:,1)-temp2r(:,2))+1,abs(temp2r(:,3)-temp2r(:,4))+1];
temp2 = temp2(:,2:3);
g = zeros(max(temp2));
x = num2cell(temp2,1);
w = sub2ind(max(temp2),x{:});
for i = 1:length(w), g(w(i)) = g(w(i)) +1; end
g = g/sum(g(:));

temp = [grpdata_syn(indgid_syn,4),inddata_syn];
temp = [temp;temp;temp;temp;temp;temp;temp;temp];
temp = temp(temp(:,1)==2,:);
% age1, age2, sex1, sex2
temp2 = [temp(2:2:end,2),temp(1:2:end,2),temp(2:2:end,3),temp(1:2:end,3)];
% [max(age), diff(age), sex1+sex2]
%temp2 = [max(temp2(:,[1,2]),[],2), abs(temp2(:,1)-temp2(:,2))+1,(temp2(:,3)-1)*2+temp2(:,4)];
temp2 = [max(temp2(:,[1,2]),[],2), abs(temp2(:,1)-temp2(:,2))+1,abs(temp2(:,3)-temp2(:,4))+1];
temp2 = temp2(:,2:3);
f = zeros(max(temp2));
x = num2cell(temp2,1);

w = sub2ind(max(temp2),x{:});

for i = 1:length(w), f(w(i)) = f(w(i)) +1; end
f = f/sum(f(:));


M = g./f;
M = max(M(:));
probm = g./(M*f);
% g is the target

% sex, sex, age difference
qq = zeros(size(temp2,1),1);
qq(:,1) = rand(size(temp2,1),1) <= probm(w);
temp1ori = [temp(2:2:end,2),temp(1:2:end,2),temp(2:2:end,3),temp(1:2:end,3)];
temp = temp(logical(reshape(repmat(qq(:,1),[1,2])',[],1)),:);
temp1r = [temp(2:2:end,2),temp(1:2:end,2),temp(2:2:end,3),temp(1:2:end,3)];


%%
close all;
h = fig('units','centimeters','width',12,'height',6,'font','helvetica','fontsize',7,'border','off');

subplot(2,1,1);
imagesc(f'/sum(f(:)));
yticks(1:2);yticklabels(0:1);
xticks(1:14);xticklabels(0:13);
xlabel('${y_1}$','interpreter','latex');
ylabel('${y_2}$','interpreter','latex')
title('\bf{Synthetic}','interpreter','latex');
caxis([0,0.4]);
subplot(2,1,2);
imagesc(g'/sum(g(:)));
yticks(1:2);yticklabels(0:1);
xticks(1:14);xticklabels(0:13);
xlabel('${y_1}$','interpreter','latex');
ylabel('${y_2}$','interpreter','latex');
title('\bf{Original (target)}','interpreter','latex');
caxis([0,0.4]);

hh = colorbar('fontsize',6);
h_bar = findobj(gcf,'Tag','Colorbar');
initpos = get(h_bar,'Position');
pos2 = [0.93,0.25,0.02,0.65];
set(h_bar, 'Position',pos2);
set(h_bar,'ticklength',0.015,'LineWidth',0.1);
set(h_bar,'ylim',[0,0.4]);

name = 'fig_rejection.pdf';
export_fig(name,'-pdf','-transparent','-painters','-rgb');
system(['start ',name]);


%%
temp00 = [grpdata(indgid,4),inddata];
temp00 = temp00(temp00(:,1)==2,:);
temp2s = [grpdata_syn(indgid_syn,4),inddata_syn];
temp2s = temp2s(temp2s(:,1)==2,:);
[y,x] = hist(temp2s(:,2),1:14); y=y/sum(y); plot(x,y,'-rs'); hold on;
[y,x] = hist(temp00(:,2),1:14); y=y/sum(y); plot(x,y,'-gs'); 

w = [temp00(1:2:end,3),temp00(2:2:end,3)];
w = [w;w(:,2),w(:,1)];
y = hist3(w,{1:2,1:2});

%%
figure;
subplot(1,2,1); imagesc(hist3(temp1ori(:,[1,2]),{1:14,1:14}));
subplot(1,2,2); imagesc(hist3(temp1r(:,[1,2]),{1:14,1:14}));


%% rejection compare
close all;
h = fig('units','centimeters','width',18,'height',5,'font','helvetica','fontsize',7,'border','off');
tempr = [grpdata(indgid,4),inddata];
tempr = tempr(tempr(:,1)==2,:);
ass1 = [tempr(1:2:end,2:end),tempr(2:2:end,2:end);
    tempr(2:2:end,2:end),tempr(1:2:end,2:end)];

ass2 = [temp(1:2:end,2:end),temp(2:2:end,2:end);
    temp(2:2:end,2:end),temp(1:2:end,2:end)];
numAtt = size(ass1,2);
cate = max(ass1);
crav = zeros(numAtt,numAtt);
for idxi = 1:numAtt-1
    for idxj = idxi+1:numAtt
        w = hist3(ass1(:,[idxi,idxj]),{1:cate(idxi),1:cate(idxj)});
        w = (w+1e-6)/sum(w(:));
        w1 = sum(w,1);
        w2 = sum(w,2);
        prd = w2*w1;
        crav(idxi,idxj) = sqrt(sum(sum((w-prd).^2./prd))/min(cate(idxi)-1,cate(idxj)-1));
    end
end
crav1 = crav+crav';
subplot(1,3,1); imagesc(crav1);  caxis([0,0.7]); hold on;
plot([5.5,5.5],[0,11],'-','color','w','linewidth',0.5); hold on;
plot([0,11],[5.5,5.5],'-','color','w','linewidth',0.5); hold on;
set(gca,'xtick',1:10,'xticklabel',{'A','S','I','L','P','A','S','I','L','P'});
set(gca,'ytick',1:10,'yticklabel',{'A','S','I','L','P','A','S','I','L','P'});
xlabel('\bf{Original}','interpreter','latex');

crav = zeros(numAtt,numAtt);
for idxi = 1:numAtt-1
    for idxj = idxi+1:numAtt
        w = hist3(ass2(:,[idxi,idxj]),{1:cate(idxi),1:cate(idxj)});
        w = (w+1e-6)/sum(w(:));
        w1 = sum(w,1);
        w2 = sum(w,2);
        prd = w2*w1;
        crav(idxi,idxj) = sqrt(sum(sum((w-prd).^2./prd))/min(cate(idxi)-1,cate(idxj)-1));
    end
end
crav2 = crav+crav';
subplot(1,3,2); imagesc(crav2); caxis([0,0.7]); hold on;
plot([5.5,5.5],[0,11],'-','color','w','linewidth',0.5); hold on;
plot([0,11],[5.5,5.5],'-','color','w','linewidth',0.5); hold on;
set(gca,'xtick',1:10,'xticklabel',{'A','S','I','L','P','A','S','I','L','P'});
set(gca,'ytick',1:10,'yticklabel',{'A','S','I','L','P','A','S','I','L','P'});
xlabel('\bf{Rejection}','interpreter','latex');

subplot(1,3,3); imagesc(abs(crav1-crav2)); caxis([0,0.7]); hold on;
plot([5.5,5.5],[0,11],'-','color','w','linewidth',0.5); hold on;
plot([0,11],[5.5,5.5],'-','color','w','linewidth',0.5); hold on;
set(gca,'xtick',1:10,'xticklabel',{'A','S','I','L','P','A','S','I','L','P'});
set(gca,'ytick',1:10,'yticklabel',{'A','S','I','L','P','A','S','I','L','P'});
xlabel('\bf{Difference}','interpreter','latex');

hh = colorbar('fontsize',6);
h_bar = findobj(gcf,'Tag','Colorbar');
initpos = get(h_bar,'Position');
pos2 = [0.93,0.25,0.02,0.65];
set(h_bar, 'Position',pos2);
set(h_bar,'ticklength',0.015,'LineWidth',0.1);
set(h_bar,'ylim',[0,0.7]);

name = 'fig_compare_two_people_rejection.pdf';
export_fig(name,'-pdf','-transparent','-painters','-rgb');
close all;
system(['start ',name]);




