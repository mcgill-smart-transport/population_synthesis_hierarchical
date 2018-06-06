% -- Lijun Sun -- %
% -- last modified: September 2, 2017 -- %

cate = indcate;
temp = num2cell(inddata_syn,1);
ind = sub2ind(cate,temp{:});
counts_ind_syn = zeros(cate);
for i = 1:size(ind,1)
    counts_ind_syn(ind(i)) = counts_ind_syn(ind(i))+1;
end

temp = num2cell(inddata,1);
ind = sub2ind(cate,temp{:});
counts_ind_raw = zeros(cate);
for i = 1:size(ind,1)
    counts_ind_raw(ind(i)) = counts_ind_raw(ind(i))+1;
end

cate = grpcate;
temp = num2cell(grpdata_syn,1);
ind = sub2ind(cate,temp{:});
counts_grp_syn = zeros(cate);
for i = 1:size(ind,1)
    counts_grp_syn(ind(i)) = counts_grp_syn(ind(i))+1;
end

temp = num2cell(grpdata,1);
ind = sub2ind(cate,temp{:});
counts_grp_raw = zeros(cate);
for i = 1:size(ind,1)
    counts_grp_raw(ind(i)) = counts_grp_raw(ind(i))+1;
end

%% marginal distribution at group level
h = fig('units','centimeters','width',18,'height',7,'font','helvetica','fontsize',7,'border','off');
labe = {'Dwell','Ethnicity','Car','Npax'};


for att = 1:length(grplvl)
    subplot(2,5,att);
    a = 1:4;
    a(att) = [];
    y1 = sumdims(counts_grp_raw,a); y1 = y1/sum(y1);
    y2 = sumdims(counts_grp_syn,a); y2 = y2/sum(y2);
    y1 = reshape(y1,[length(y1),1]);
    y2 = reshape(y2,[length(y2),1]);
    b = bar(1:grpcate(att),[y1,y2]);
    b(1).FaceColor = scicol('red');
    b(1).EdgeColor = 'none';
    b(2).FaceColor = scicol('green');
    b(2).EdgeColor = 'none';
    xlim([0.5,grpcate(att)+0.5]);
    xlabel(labe{att});
    if att == 4
        legendflex(b, {'Original','Synthetic'}, 'ref', gca, 'nrow',2, 'anchor', {'ne','ne'},'buffer',[80 -5],'box','off','xscale',0.2);
    end
    set(gca,'Box','on','TickDir', 'in','TickLength',[.02 .05],'layer','top','color','none',...
        'XMinorTick','off','YMinorTick','on','ZMinorTick','on','YGrid','off','XColor',[.15,.15,.15],'YColor',[.15,.15,.15],'ZColor',[.15,.15,.15],'LineWidth',0.5);
end

labe = {'Age','Sex','Income','License','Pass'};
for att = 1:length(indlvl)
    subplot(2,5,att+5);
    a = 1:5;
    a(att) = [];
    y1 = sumdims(counts_ind_raw,a); y1 = y1/sum(y1);
    y2 = sumdims(counts_ind_syn,a); y2 = y2/sum(y2);
    y1 = reshape(y1,[length(y1),1]);
    y2 = reshape(y2,[length(y2),1]);
    
    b = bar(1:indcate(att),[y1,y2]);
    b(1).FaceColor = scicol('red');
    b(1).EdgeColor = 'none';
    b(2).FaceColor = scicol('green');
    b(2).EdgeColor = 'none';
    xlim([0.5,indcate(att)+0.5]);
    xlabel(labe{att});
    if indcate(att)>10
        set(gca,'xtick',1:2:indcate(att));
    end
    set(gca,'Box','on','TickDir', 'in','TickLength',[.02 .05],'layer','top','color','none',...
        'XMinorTick','off','YMinorTick','on','ZMinorTick','on','YGrid','off','XColor',[.15,.15,.15],'YColor',[.15,.15,.15],'ZColor',[.15,.15,.15],'LineWidth',0.5);
end

name = 'fig_marginal.eps';
export_fig(name,'-pdf','-transparent','-painters','-rgb');
close all;
system(['start ',name]);

%%
clc;
h = fig('units','centimeters','width',18,'height',7,'font','helvetica','fontsize',7,'border','off');
labe = {'Dwell','Ethnicity','Car','Npax'};



subplot(2,5,1);
t1 = counts_grp_raw;
t2 = counts_grp_syn;
t1 = t1/sum(t1(:));
t2 = t2/sum(t2(:));

m = 0;
plot([0,1],[0,1],'color',scicol('green'),'linewidth',1); hold on;
for att = 1:length(grpcate)
    a = 1:4;
    a(att) = [];
    d1 = sumdims(t1,a);
    d2 = sumdims(t2,a);
    d1 = d1(:);
    d2 = d2(:);
    plot(d1(d1>1e-3),d2(d1>1e-3),'o','color',scicol('red'),'markersize',2); hold on;
    m = max(m,max(d1));
end
xlabel('\bf{Original}','interpreter','latex');
ylabel('$\bf{Synthetic}$ (household)','interpreter','latex');
text(0.1*m,m,'1D');
xlim([0,m*1.1]);
ylim([0,m*1.1]);


m = 0;
subplot(2,5,2);
plot([0,1],[0,1],'color',scicol('green'),'linewidth',1); hold on;
for att1 = 1:length(grpcate)
    for att2 = 1:length(grpcate)
        a = 1:4;
        if length(unique([att1,att2])) == 2
            a([att1,att2]) = [];
            d1 = sumdims(t1,a);
            d2 = sumdims(t2,a);
            d1 = d1(:);
            d2 = d2(:);
            plot(d1(d1>1e-3),d2(d1>1e-3),'o','color',scicol('red'),'markersize',2); hold on;
            m = max(m,max(d1(:)));
        end
    end
end

xlabel('\bf{Original}','interpreter','latex');
text(0.1*m,m,'2D');
xlim([0,m*1.1]);
ylim([0,m*1.1]);



m = 0;
subplot(2,5,3);
plot([0,1],[0,1],'color',scicol('green'),'linewidth',1); hold on;
for att1 = 1:length(grpcate)
    for att2 = 1:length(grpcate)
        for att3 = 1:length(grpcate)
            if length(unique([att1,att2,att3])) == 3
                a = 1:4;
                a([att1,att2,att3]) = [];
                d1 = sumdims(t1,a);
                d2 = sumdims(t2,a);
                d1 = d1(:);
                d2 = d2(:);
                plot(d1(d1>1e-3),d2(d1>1e-3),'o','color',scicol('red'),'markersize',2); hold on;
                m = max(m,max(d1(:)));
            end
        end
    end
end

xlabel('\bf{Original}','interpreter','latex');
text(0.1*m,m,'3D');
xlim([0,m*1.1]);
ylim([0,m*1.1]);


m = 0;
subplot(2,5,4);
plot([0,1],[0,1],'color',scicol('green'),'linewidth',1); hold on;

d1 = t1;
d2 = t2;
d1 = d1(:);
d2 = d2(:);
plot(d1(d1>1e-3),d2(d1>1e-3),'o','color',scicol('red'),'markersize',2); hold on;
m = max(m,max(d1(:)));


xlabel('\bf{Original}','interpreter','latex');
text(0.1*m,m,'4D');
xlim([0,m*1.1]);
ylim([0,m*1.1]);



subplot(2,5,6);
t1 = counts_ind_raw;
t2 = counts_ind_syn;
t1 = t1/sum(t1(:));
t2 = t2/sum(t2(:));

m = 0;
plot([0,1],[0,1],'color',scicol('green'),'linewidth',1); hold on;
for att = 1:length(indcate)
    a = 1:5;
    a(att) = [];
    d1 = sumdims(t1,a);
    d2 = sumdims(t2,a);
    d1 = d1(:);
    d2 = d2(:);
    plot(d1(d1>1e-3),d2(d1>1e-3),'o','color',scicol('red'),'markersize',2); hold on;
    m = max(m,max(d1));
end
xlabel('\bf{Original}','interpreter','latex');
ylabel('$\bf{Synthetic}$ (individual)','interpreter','latex');
text(0.1*m,m,'1D');
xlim([0,m*1.1]);
ylim([0,m*1.1]);


m = 0;
subplot(2,5,7);
plot([0,1],[0,1],'color',scicol('green'),'linewidth',1); hold on;
for att1 = 1:length(indcate)
    for att2 = 1:length(indcate)
        a = 1:5;
        if length(unique([att1,att2])) == 2
            a([att1,att2]) = [];
            d1 = sumdims(t1,a);
            d2 = sumdims(t2,a);
            d1 = d1(:);
            d2 = d2(:);
            plot(d1(d1>1e-3),d2(d1>1e-3),'o','color',scicol('red'),'markersize',2); hold on;
            m = max(m,max(d1(:)));
        end
    end
end

xlabel('\bf{Original}','interpreter','latex');
text(0.1*m,m,'2D');
xlim([0,m*1.1]);
ylim([0,m*1.1]);



m = 0;
subplot(2,5,8);
plot([0,1],[0,1],'color',scicol('green'),'linewidth',1); hold on;
for att1 = 1:length(indcate)
    for att2 = 1:length(indcate)
        for att3 = 1:length(indcate)
            if length(unique([att1,att2,att3])) == 3
                a = 1:5;
                a([att1,att2,att3]) = [];
                d1 = sumdims(t1,a);
                d2 = sumdims(t2,a);
                d1 = d1(:);
                d2 = d2(:);
                plot(d1(d1>1e-3),d2(d1>1e-3),'o','color',scicol('red'),'markersize',2); hold on;
                m = max(m,max(d1(:)));
            end
        end
    end
end

xlabel('\bf{Original}','interpreter','latex');
text(0.1*m,m,'3D');
xlim([0,m*1.1]);
ylim([0,m*1.1]);


m = 0;
subplot(2,5,9);
plot([0,1],[0,1],'color',scicol('green'),'linewidth',1); hold on;
for att1 = 1:length(indcate)
    for att2 = 1:length(indcate)
        for att3 = 1:length(indcate)
            for att4 = 1:length(indcate)
                if length(unique([att1,att2,att3,att4])) == 4
                    a = 1:5;
                    a([att1,att2,att3,att4]) = [];
                    d1 = sumdims(t1,a);
                    d2 = sumdims(t2,a);
                    d1 = d1(:);
                    d2 = d2(:);
                    plot(d1(d1>1e-3),d2(d1>1e-3),'o','color',scicol('red'),'markersize',2); hold on;
                    m = max(m,max(d1(:)));
                end
            end
        end
    end
end

xlabel('\bf{Original}','interpreter','latex');
text(0.1*m,m,'4D');
xlim([0,m*1.1]);
ylim([0,m*1.1]);


m = 0;
subplot(2,5,10);
plot([0,1],[0,1],'color',scicol('green'),'linewidth',1); hold on;

d1 = t1;
d2 = t2;
d1 = d1(:);
d2 = d2(:);
plot(d1(d1>1e-3),d2(d1>1e-3),'o','color',scicol('red'),'markersize',2); hold on;
m = max(m,max(d1(:)));


xlabel('\bf{Original}','interpreter','latex');
text(0.1*m,m,'5D');
xlim([0,m*1.1]);
ylim([0,m*1.1]);



name = 'fig_tucker_individual.eps';
export_fig(name,'-pdf','-transparent','-painters','-rgb');
close all;
system(['start ',name]);


%% joint association at both group and individual levels
close all;
h = fig('units','centimeters','width',12,'height',5,'font','helvetica','fontsize',6,'border','off');
temp = [grpdata(indgid,:),inddata];
cate = max(temp,[],1);
numAtt = size(temp,2);
crav = zeros(numAtt,numAtt);
for idxi = 1:numAtt-1
    for idxj = idxi+1:numAtt
        w = hist3(temp(:,[idxi,idxj]),{1:cate(idxi),1:cate(idxj)});
        w = (w+1e-6)/sum(sum(w));
        w1 = sum(w,1);
        w2 = sum(w,2);
        prd = w2*w1;
        crav(idxi,idxj) = sqrt(sum(sum((w-prd).^2./prd))/min(cate(idxi)-1,cate(idxj)-1));
    end
end
crav1 = crav + crav';
disp(max(max(crav1)));
subplot(1,2,1);
imagesc(crav1); caxis([0,0.5]);hold on;
plot([4.5,4.5],[0,11],'-','color','w','linewidth',1); hold on;
plot([0,11],[4.5,4.5],'-','color','w','linewidth',1); hold on;

set(gca,'xtick',1:numAtt,'xticklabel',{'Dwell','Eth','Car','Npax','Age','Sex','Inc','Lic','Pass'});
set(gca,'ytick',1:numAtt,'yticklabel',{'Dwell','Eth','Car','Npax','Age','Sex','Inc','Lic','Pass'});
xlabel('\bf{Original}','interpreter','latex');
%text(-0.25,1.05,'(a)','unit','normalized','VerticalAlignment','top','fontname','times','fontsize',10,'fontweight','bold');



temp = [grpdata_syn(indgid_syn,:),inddata_syn];
crav = zeros(numAtt,numAtt);
for idxi = 1:numAtt-1
    for idxj = idxi+1:numAtt
        w = hist3(temp(:,[idxi,idxj]),{1:cate(idxi),1:cate(idxj)});
        w = (w+1e-6)/sum(sum(w));
        w1 = sum(w,1);
        w2 = sum(w,2);
        prd = w2*w1;
        crav(idxi,idxj) = sqrt(sum(sum((w-prd).^2./prd))/min(cate(idxi)-1,cate(idxj)-1));
    end
end
crav2 = crav + crav';
disp(max(max(crav2)));
subplot(1,2,2);
imagesc(crav2); caxis([0,0.5]); hold on;
plot([4.5,4.5],[0,11],'-','color','w','linewidth',1); hold on;
plot([0,11],[4.5,4.5],'-','color','w','linewidth',1); hold on;

set(gca,'xtick',1:numAtt,'xticklabel',{'Dwell','Eth','Car','Npax','Age','Sex','Inc','Lic','Pass'});
set(gca,'ytick',1:numAtt,'yticklabel',{'Dwell','Eth','Car','Npax','Age','Sex','Inc','Lic','Pass'});
xlabel('\bf{Synthetic}','interpreter','latex');
%text(-0.25,1.05,'(b)','unit','normalized','VerticalAlignment','top','fontname','times','fontsize',10,'fontweight','bold');



hh = colorbar('fontsize',6);
h_bar = findobj(gcf,'Tag','Colorbar');
initpos = get(h_bar,'Position');
pos2 = [0.95,0.15,0.02,0.7];
set(h_bar, 'Position',pos2);
set(h_bar,'ticklength',0.015,'LineWidth',0.1);
set(h_bar,'ylim',[0,0.5]);


name = 'fig_compare_joint.eps';
export_fig(name,'-pdf','-transparent','-painters','-cmyk');
%close all;
system(['start ',name]);


%%
w1 = crav1(1:4,1:4);
w2 = crav2(1:4,1:4);
plot(w1(:),w2(:),'s')


figure;
w1 = crav1(5:9,5:9);
w2 = crav2(5:9,5:9);
plot(w1(:),w2(:),'s')



%%
