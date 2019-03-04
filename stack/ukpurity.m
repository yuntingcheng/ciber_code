flight=40030;
inst=1;
ifield=8;
alpha=-15;
beta=258;
m_max=17;
catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/UKIDSS/');
savedir='/Users/ytcheng/ciber/doc/20170617_Stacking/plots/ukpurity/';

quad_arr=['A','B','C','D'];
dt=get_dark_times(flight,inst,ifield);

strmask=zeros(1024);

iquad=1;
%%% read cat data %%%
quad=quad_arr(iquad);
catfile=strcat(catdir,dt.name,'_',quad,'_uk.txt');
M = csvread(catfile,1);
x_arr=squeeze(M(:,5)');
y_arr=squeeze(M(:,4)');
if inst==1
    m_arr=squeeze(M(:,6)');
else
    m_arr=squeeze(M(:,7)');
end

if ifield==8
    cls_arr=squeeze(M(:,11)');
    ps_arr=squeeze(M(:,12)');
    pg_arr=squeeze(M(:,13)');
    pn_arr=squeeze(M(:,14)');
    psat_arr=squeeze(M(:,15)');
else
    cls_arr=squeeze(M(:,12)');
    ps_arr=squeeze(M(:,13)');
    pg_arr=squeeze(M(:,14)');
    pn_arr=squeeze(M(:,15)');
    psat_arr=squeeze(M(:,16)');    
end
%%
figure
binedge=0:5e-3:1;
bin=(binedge(2:end)+binedge(1:end-1))/2;
[hval,edges] = histcounts(pn_arr,binedge);
sp=find(hval~=0);
hval=hval(sp);bin=bin(sp);
bar(bin,log10(hval));hold on
plot(bin,log10(hval),'o');
xlim([0,0.6])
line([0.1 0.1], get(gca, 'ylim'),'color','r');
line([0.06 0.06], get(gca, 'ylim'),'color','r');
line([0.04 0.04], get(gca, 'ylim'),'color','r');
line([0.02 0.02], get(gca, 'ylim'),'color','r');
title(strcat(dt.name,'\_',quad))
xlabel('pNoise','fontsize',20)
ylabel('log10(counts)','fontsize',20)
savename=strcat(savedir,dt.name,'_pnhist');
print(savename,'-dpng');

%%
ukmapt=zeros(1024);
ukmaps=zeros(1024);
ukmapg=zeros(1024);
ukmap1=zeros(1024);
ukmap2=zeros(1024);

for iquad=1:4
quad=quad_arr(iquad);
srcmapdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
sukmapt=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapt.fits'));
sukmaps=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmaps.fits'));
sukmapg=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapg.fits'));
sukmap1=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmap1.fits'));
sukmap2=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmap2.fits'));

if iquad==1
    ukmapt(1:512,1:512)=sukmapt;
    ukmaps(1:512,1:512)=sukmaps;
    ukmapg(1:512,1:512)=sukmapg;
    ukmap1(1:512,1:512)=sukmap1;
    ukmap2(1:512,1:512)=sukmap2;
elseif iquad==2
    ukmapt(513:1024,1:512)=sukmapt;
    ukmaps(513:1024,1:512)=sukmaps;
    ukmapg(513:1024,1:512)=sukmapg;
    ukmap1(513:1024,1:512)=sukmap1;
    ukmap2(513:1024,1:512)=sukmap2;
elseif iquad==3
    ukmapt(1:512,513:1024)=sukmapt;
    ukmaps(1:512,513:1024)=sukmaps;
    ukmapg(1:512,513:1024)=sukmapg;
    ukmap1(1:512,513:1024)=sukmap1;
    ukmap2(1:512,513:1024)=sukmap2;
elseif iquad==4
    ukmapt(513:1024,513:1024)=sukmapt;
    ukmaps(513:1024,513:1024)=sukmaps;
    ukmapg(513:1024,513:1024)=sukmapg;
    ukmap1(513:1024,513:1024)=sukmap1;
    ukmap2(513:1024,513:1024)=sukmap2;
end
end

%%
figure
setwinsize(gcf,1000,800)
imageclip(log10(ukmapt));
title(strcat(dt.name,' (stars+galaxies+noise)'),'fontsize',20);
h = colorbar;
ylabel(h, 'log10($\nu I_\nu)[nW/m^2/sr]$','interpreter','latex','fontsize',20);
v=caxis;
savename=strcat(savedir,dt.name,'_map_all');
print(savename,'-dpng');
%%
figure
setwinsize(gcf,1000,800)


subplot(2,2,1)
imageclip(log10(ukmaps));
title('stars (pNoise<0.04, pStar>pGalaxy)','fontsize',20);
h = colorbar;
ylabel(h, 'log10($\nu I_\nu)[nW/m^2/sr]$','interpreter','latex','fontsize',20);
caxis(v);

caxis(v);
subplot(2,2,2)
imageclip(log10(ukmapg));
title('galaxies (pNoise<0.04, pStar<pGalaxy)','fontsize',20);
h = colorbar;
ylabel(h, 'log10($\nu I_\nu)[nW/m^2/sr]$','interpreter','latex','fontsize',20);
caxis(v);

subplot(2,2,3)
imageclip(log10(ukmap1));
title('pNoise>0.1','fontsize',20);
h = colorbar;
ylabel(h, 'log10($\nu I_\nu)[nW/m^2/sr]$','interpreter','latex','fontsize',20);
caxis(v);

subplot(2,2,4)
imageclip(log10(ukmap2));
title('0.04<pNoise<0.1','fontsize',20);
h = colorbar;
ylabel(h, 'log10($\nu I_\nu)[nW/m^2/sr]$','interpreter','latex','fontsize',20);
caxis(v);
savename=strcat(savedir,dt.name,'_map_components');
print(savename,'-dpng');

%%
figure
setwinsize(gcf,1000,800)
imageclip(log10(ukmaps+ukmapg));
title(strcat(dt.name,' (stars+galaxies)'),'fontsize',20);
h = colorbar;
ylabel(h, 'log10($\nu I_\nu)[nW/m^2/sr]$','interpreter','latex','fontsize',20);
caxis(v);
savename=strcat(savedir,dt.name,'_map_sg');
print(savename,'-dpng');

%%
if ifield~=8
    [~,masks]=get_skymap(ukmaps,ones(1024),5,2);
    [~,maskg]=get_skymap(ukmapg,ones(1024),5,2);
else
    strmask=ones(1024);
    strmask(1:120,800:end)=0;
    strmask(980:end,800:end)=0;
    [~,masks]=get_skymap(ukmaps,strmask,5,2);
    [~,maskg]=get_skymap(ukmapg,strmask,5,2); 
end
%%

[sms]=fillpadsmooth(ukmaps,masks,30);
[smg]=fillpadsmooth(ukmapg,maskg,30);

figure
setwinsize(gcf,1000,800)
subplot(2,2,1)
imageclip(log10(ukmaps));
title('Stars','fontsize',20);
h = colorbar;
ylabel(h, 'log10($\nu I_\nu)[nW/m^2/sr]$','interpreter','latex','fontsize',20);
subplot(2,2,2)
imageclip(log10(ukmapg));
title('Galaxies','fontsize',20);
h = colorbar;
ylabel(h, 'log10($\nu I_\nu)[nW/m^2/sr]$','interpreter','latex','fontsize',20);
subplot(2,2,3)
imageclip(sms);
title('Stars smoothed','fontsize',20);
h = colorbar;
ylabel(h, '$\nu I_\nu[nW/m^2/sr]$','interpreter','latex','fontsize',20);
subplot(2,2,4)
imageclip(smg);
title('Galaxeis smoothed','fontsize',20);
h = colorbar;
ylabel(h, '$\nu I_\nu[nW/m^2/sr]$','interpreter','latex','fontsize',20);

savename=strcat(savedir,dt.name,'_map_sm');
print(savename,'-dpng');

%%
mask=masks.*maskg;

ukmaps1=ukmaps.*mask;
ukmaps1=ukmaps1-mean(ukmaps1(find(mask)));
ukmaps1=ukmaps1.*mask;

ukmapg1=ukmapg.*mask;
ukmapg1=ukmapg1-mean(ukmapg1(find(mask)));
ukmapg1=ukmapg1.*mask;


[Cls,l]=get_angular_spec(ukmaps1,ukmaps1,7);
[Clg,l]=get_angular_spec(ukmapg1,ukmapg1,7);
[Clx,l]=get_angular_spec(ukmaps1,ukmapg1,7);

loglog(l,l.*(l+1).*Cls./2./pi,'Displayname','star');hold on
loglog(l,l.*(l+1).*Clg./2./pi,'Displayname','gal');hold on
loglog(l,l.*(l+1).*Clx./2./pi,'Displayname','star x gal');hold on
xlim([1e2,2e5]);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
leg=legend('show','location','northwest');
legend boxoff
set(leg,'FontSize',15);
title(dt.name,'FontSize',15);
savename=strcat(savedir,dt.name,'_Cl');
print(savename,'-dpng');
%%
xr=Clx./sqrt(Clg.*Cls);

semilogx(l,xr,'linewidth',2);hold on
plot([1e2,2e5],[0,0],'k');
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$r=C\ell_x/\sqrt{C\ell_sC\ell_g}$',...
        'interpreter','latex','fontsize',18)
title(dt.name,'FontSize',15);
xlim([1e2,2e5]);
ylim([-0.4,0.4]);
savename=strcat(savedir,dt.name,'_xcoeff');
print(savename,'-dpng');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot UK SWIRE different m bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
ifield=8;
dt=get_dark_times(flight,inst,ifield);

ukmapt=zeros(1024);
ukmaps=zeros(1024);
ukmapg=zeros(1024);
ukmap1=zeros(1024);
ukmap2=zeros(1024);
quad_arr=['A','B','C','D'];

for iquad=1:4
quad=quad_arr(iquad);
srcmapdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
sukmapt=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapt.fits'));
sukmaps=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmaps.fits'));
sukmapg=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapg.fits'));
sukmap1=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmap1.fits'));
sukmap2=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmap2.fits'));

% srcmapt some how is wrong
sukmapt=sukmaps+sukmapg+sukmap1+sukmap2;
fits_write(strcat(srcmapdir,dt.name,'_',quad,'_srcmapt.fits'),sukmapt);

if iquad==1
    ukmapt(1:512,1:512)=sukmapt;
    ukmaps(1:512,1:512)=sukmaps;
    ukmapg(1:512,1:512)=sukmapg;
    ukmap1(1:512,1:512)=sukmap1;
    ukmap2(1:512,1:512)=sukmap2;
elseif iquad==2
    ukmapt(513:1024,1:512)=sukmapt;
    ukmaps(513:1024,1:512)=sukmaps;
    ukmapg(513:1024,1:512)=sukmapg;
    ukmap1(513:1024,1:512)=sukmap1;
    ukmap2(513:1024,1:512)=sukmap2;
elseif iquad==3
    ukmapt(1:512,513:1024)=sukmapt;
    ukmaps(1:512,513:1024)=sukmaps;
    ukmapg(1:512,513:1024)=sukmapg;
    ukmap1(1:512,513:1024)=sukmap1;
    ukmap2(1:512,513:1024)=sukmap2;
elseif iquad==4
    ukmapt(513:1024,513:1024)=sukmapt;
    ukmaps(513:1024,513:1024)=sukmaps;
    ukmapg(513:1024,513:1024)=sukmapg;
    ukmap1(513:1024,513:1024)=sukmap1;
    ukmap2(513:1024,513:1024)=sukmap2;
end
end
%%
m_arr=17:24;
setwinsize(gcf,6000,750)

subplot(3,3,1)
imageclip(log10(ukmapg));
title('total galaxies','fontsize',20)
h = colorbar;
ylabel(h, '$log_{10}(\nu I_\nu)[nW/m^2/sr]$','interpreter','latex','fontsize',15);
v=caxis;
for m=m_arr
    ukmap=zeros(1024);
    for iquad=1:4
            quad=quad_arr(iquad);
            sukmap=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapg',...
                num2str(m),'.fits'));
        if iquad==1
            ukmap(1:512,1:512)=sukmap;
        elseif iquad==2
            ukmap(513:1024,1:512)=sukmap;
        elseif iquad==3
            ukmap(1:512,513:1024)=sukmap;
        elseif iquad==4
            ukmap(513:1024,513:1024)=sukmap;
        end

    end
    
    
    subplot(3,3,m-15)
    imageclip(log10(ukmap));
    h = colorbar;
    ylabel(h, '$log_{10}(\nu I_\nu)[nW/m^2/sr]$',...
        'interpreter','latex','fontsize',15);

    if m==17
        title('m<17','fontsize',20);
    elseif m==24
        title('m>23','fontsize',20);
    else
        title(strcat(num2str(m),'<m<',num2str(m+1)),'fontsize',20);
    end
    
    caxis(v)
    
end
savename=strcat(savedir,dt.name,'_mag_g');
print(savename,'-dpng');
