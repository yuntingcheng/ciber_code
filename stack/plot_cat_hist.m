%% compare 2M and UK
flight=40030;
inst=2;
ifield=5;
catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/');
quad='C';
dt=get_dark_times(flight,inst,ifield);

catfile1=strcat(catdir,'PSC/',dt.name,'_',quad,'_2m.txt');
catfile2=strcat(catdir,'UKIDSS/',dt.name,'_',quad,'_uk.txt');

M2m = csvread(catfile1,1);
Muk = csvread(catfile2,1);
if inst==1
    m2m_arr=squeeze(M2m(:,6)');
    muk_arr=squeeze(Muk(:,6)');
else
    m2m_arr=squeeze(M2m(:,7)');
    muk_arr=squeeze(Muk(:,7)');
end


binedge=1:0.2:24;
histogram(muk_arr,'BinEdges',binedge);hold on
histogram(m2m_arr,'BinEdges',binedge);
legend({'UKIDSS','2MASS'},'location','northwest')
xlabel('AB mag (CIBER I band)')
ylabel('counts')
title(strcat(dt.name,'\_',quad))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare different cats in 2M
flight=40030;
inst=1;
ifield=4;
dt=get_dark_times(flight,inst,ifield);
quad='A';
catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/');

binedge=6:0.25:24;
bin=(binedge(2:end)+binedge(1:end-1))./2;

catname='PSC';
catfile=strcat(catdir,catname,'/',dt.name,'_',quad,'_2m.txt');
M = csvread(catfile,1);
if inst==1
    m_arr=squeeze(M(:,6)');
else
    m_arr=squeeze(M(:,7)');
end
psccounts=numel(m_arr);
hpsc=histcounts(m_arr,'BinEdges',binedge,'Normalization','Probability');

catname='XSC';
catfile=strcat(catdir,catname,'/',dt.name,'_',quad,'_2m.txt');
M = csvread(catfile,1);
if inst==1
    m_arr=squeeze(M(:,6)');
else
    m_arr=squeeze(M(:,7)');
end
xsccounts=numel(m_arr);
hxsc=histcounts(m_arr,'BinEdges',binedge,'Normalization','Probability');

catname='XSCrej';
catfile=strcat(catdir,catname,'/',dt.name,'_',quad,'_2m.txt');
M = csvread(catfile,1);
if inst==1
    m_arr=squeeze(M(:,6)');
else
    m_arr=squeeze(M(:,7)');
end
xscrejcounts=numel(m_arr);
hxscrej=histcounts(m_arr,'BinEdges',binedge,'Normalization','Probability');

rel_arr=['A','B','C','D','E','F'];
for irel=1:6
    rel=rel_arr(irel);
    catname='PSCrej';
    catfile=strcat(catdir,catname,'/',dt.name,'_',quad,'_2m_',rel,'.txt');
    M = csvread(catfile,1);
    if inst==1
        m_arr=squeeze(M(:,6)');
    else
        m_arr=squeeze(M(:,7)');
    end
    h=histcounts(m_arr,'BinEdges',binedge,'Normalization','Probability');
    pscrej(irel).counts=numel(m_arr);
    pscrej(irel).h=h;
end

savedir='/Users/ytcheng/ciber/doc/20170617_Stacking/plots/cats/';
counts_arr=[psccounts,xsccounts,pscrej(1).counts,pscrej(2).counts,...
 pscrej(3).counts,pscrej(4).counts,pscrej(5).counts,pscrej(6).counts,xscrejcounts];
Xt=1:length(counts_arr);
semilogy(Xt,counts_arr,'.','markersize',20);


set(gca,'XTick',Xt);
algos = ['    PSC'; '    XSC'; 'PSCrejA'; 'PSCrejB'; 'PSCrejC'; 'PSCrejD';
     'PSCrejE'; 'PSCrejF'; ' XSCrej'];

ax = axis; % Current axis limits
axis(axis); % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4); % Y-axis limits

% Remove the default labels
set(gca,'XTickLabel','')

% Place the text labels
t = text(Xt,Yl(1)*ones(1,length(Xt)),algos(1:length(Xt),:));
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
'Rotation',45, 'Fontsize', 13);
ylabel('counts', 'Fontsize', 15)

xlim([0.5,9.5])
title('2MASS catalog counts (elat10\_A)')

savename=strcat(savedir,'counts');
print(savename,'-dpng');

plot(bin,hpsc.*4,'k','linewidth',2);hold on
plot(bin,hxsc.*4,'b','linewidth',2);hold on
xlabel('AB mag','fontsize',15);
ylabel('frac coutns/AB mag','fontsize',15);
legend({'PSC','XSC'},'fontsize',15);
legend boxoff
title('2MASS catalog CIBER I-band mag hist (elat10\_A)')
savename=strcat(savedir,'hist_xsc');
print(savename,'-dpng');

plot(bin,hpsc.*4,'k','linewidth',2);hold on
plot(bin,pscrej(1).h.*4,'b','linewidth',2);hold on
plot(bin,pscrej(6).h.*4,'r','linewidth',2);hold on
xlabel('AB mag','fontsize',15);
ylabel('frac coutns/AB mag','fontsize',15);
legend({'PSC','PSC rej A','PSC rej F'},'location','northwest','fontsize',15);
legend boxoff
title('2MASS catalog CIBER I-band mag hist (elat10\_A)')
savename=strcat(savedir,'hist_pscrej');
print(savename,'-dpng');
%%  UK star and gal compare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savedir='/Users/ytcheng/ciber/doc/20170617_Stacking/plots/cats/';
flight=40030;
inst=2;
ifield=8;
dt=get_dark_times(flight,inst,ifield);
quad='A';
catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/UKIDSS/');
catfile=strcat(catdir,dt.name,'_',quad,'_uk.txt');
M = csvread(catfile,1);
if inst==1
    m_arr=squeeze(M(:,6)');
else
    m_arr=squeeze(M(:,7)');
end

ps_arr=squeeze(M(:,12)');
pg_arr=squeeze(M(:,13)');
pn_arr=squeeze(M(:,14)');

sp_s=find(pn_arr<0.04 & ps_arr>pg_arr);
sp_g=find(pn_arr<0.04 & pg_arr>ps_arr);

dm=0.2;
binedge=9.5:dm:25.5;
mbin_arr=(binedge(2:end)+binedge(1:end-1))/2;
Ns_arr = histcounts(m_arr(sp_s),binedge);Ns_arr=Ns_arr./dm;
Ng_arr = histcounts(m_arr(sp_g),binedge);Ng_arr=Ng_arr./dm;
plot(mbin_arr,Ns_arr,'r','linewidth',2);hold on
plot(mbin_arr,Ng_arr,'b','linewidth',2);
xlabel('AB mag','fontsize',15);
ylabel('dN/dm','fontsize',15);
legend({'stars','galaxies'},'location','northwest','fontsize',15);
legend boxoff
title(strcat(dt.name,'\_',quad,' band',num2str(inst)))
savename=strcat(savedir,'hist_uksg_TM',num2str(inst));
print(savename,'-dpng');
%% compare UK,2m,trilegal
% UKIDSS
%pncut=0.04;
pncut=0.0031;
figure
savedir='/Users/ytcheng/ciber/doc/20170617_Stacking/plots/cats/';
flight=40030;
inst=1;
ifield=8;
dt=get_dark_times(flight,inst,ifield);
quad='A';
catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/UKIDSS/');
catfile=strcat(catdir,dt.name,'_',quad,'_uk.txt');
M = csvread(catfile,1);
x_arr=squeeze(M(:,5)');
y_arr=squeeze(M(:,4)');
x_arr=x_arr+1;
y_arr=y_arr+1;

if inst==1
    m_arr=squeeze(M(:,6)');
else
    m_arr=squeeze(M(:,7)');
end
cls_arr=squeeze(M(:,11)');
ps_arr=squeeze(M(:,12)');
pg_arr=squeeze(M(:,13)');
pn_arr=squeeze(M(:,14)');

sp_s=find(pn_arr<pncut & cls_arr==-1 ...
    & x_arr>0.5 & x_arr<512.5 & y_arr>0.5 & y_arr<512.5);
sp_g=find(pn_arr<pncut & cls_arr==1 ...
    & x_arr>0.5 & x_arr<512.5 & y_arr>0.5 & y_arr<512.5);

sp = find(x_arr>0.5 & x_arr<512.5 & y_arr>0.5 & y_arr<512.5);
dm=0.4;
binedge=9.5:dm:25.5;
mbin_arr=(binedge(2:end)+binedge(1:end-1))/2;
Ns_arr = histcounts(m_arr(sp_s),binedge);Ns_arr=Ns_arr./dm;
Ng_arr = histcounts(m_arr(sp_g),binedge);Ng_arr=Ng_arr./dm;
N_arr = histcounts(m_arr(sp),binedge);N_arr=N_arr./dm;
plot(mbin_arr,Ns_arr,'r','linewidth',2);hold on
plot(mbin_arr,Ng_arr,'b','linewidth',2);
%plot(mbin_arr,N_arr,'m','linewidth',2);

% 2MASS PSC
catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/');
catname='PSC';
catfile=strcat(catdir,catname,'/',dt.name,'_',quad,'_2m.txt');
M = csvread(catfile,1);
x_arr=squeeze(M(:,5)');
y_arr=squeeze(M(:,4)');

x_arr=x_arr+1;
y_arr=y_arr+1;

if inst==1
    m_arr=squeeze(M(:,6)');
else
    m_arr=squeeze(M(:,7)');
end

sp=find(x_arr>0.5 & x_arr<512.5 & y_arr>0.5 & y_arr<512.5);
y_arr=squeeze(M(:,4)');

x_arr=x_arr+1;
y_arr=y_arr+1;

N2m_arr=histcounts(m_arr(sp),binedge);N2m_arr=N2m_arr./dm;
plot(mbin_arr,N2m_arr,'g','linewidth',2);

% Trilegal
catdir='/Users/ytcheng/ciber/doc/20171018_stackihl/trilegal/';
catfile=strcat(catdir,dt.name,'_',quad,'_srcmapcoord.txt');
M = csvread(catfile);
m_arr=M(:,1);
Ntri_arr=histcounts(m_arr,binedge);Ntri_arr=Ntri_arr./dm;
plot(mbin_arr,Ntri_arr,'k','linewidth',2);
plot(mbin_arr,Ntri_arr.*1.25,'k--','linewidth',2);

set(gca,'YScale', 'log');
xlabel('AB mag','fontsize',15);
ylabel('dN/dm','fontsize',15);
legend({'UKstars','UKgalaxies','2mass','trilegal','trilegal x1.25'},...
    'location','northwest','fontsize',15);
legend boxoff
title(strcat(dt.name,'\_',quad,' band',num2str(inst),' (lowest pNoise)'))
savename=strcat(savedir,'hist_trilegal_TM',num2str(inst));
print(savename,'-dpng');
