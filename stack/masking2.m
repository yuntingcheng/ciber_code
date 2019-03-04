flight=40030;
mypaths=get_paths(flight);
inst=1;
pixscale=7;

cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;

quad_arr=['A','B','C','D'];
srcmapdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
load(strcat(mypaths.alldat,'TM',num2str(inst),'/','maskdat'),'maskdat');
load(strcat(mypaths.alldat,'TM',num2str(inst),'/','FFdat'),'FFdat');

load(sprintf('%sTM%d/darklongdat',mypaths.filtmap,inst),'darklongdat');
DCtemplate=darklongdat.DCtemplate; clear darklongdat

savedir='/Users/ytcheng/ciber/doc/20170617_Stacking/plots/cats/';
%%
ifield=4;
dt=get_dark_times(flight,inst,ifield);

FF=FFdat(ifield).FF;
%%% get flight map %%%
loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
load(strcat(loaddir,'flightmap'),'flightmap');
rawmap=flightmap.filtmapf;
calmap=(rawmap-DCtemplate)./FF;
calmap(find(calmap~=calmap))=0;
calmap(find(calmap==-inf))=0;
calmap(find(calmap==inf))=0;
%%
%%% get sim srcmap %%%
ukmap=zeros(1024);
tmmap=zeros(1024);
for iquad=1:4
    quad=quad_arr(iquad);
    sukmap=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapt.fits'));
    %sukmaps=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmaps.fits'));
    %sukmapg=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapg.fits'));
    %sukmap=sukmaps+sukmapg;
    stmmap=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapt_2m.fits'));
    
    if iquad==1
        ukmap(1:512,1:512)=sukmap;
        tmmap(1:512,1:512)=stmmap;
    elseif iquad==2
        ukmap(513:1024,1:512)=sukmap;
        tmmap(513:1024,1:512)=stmmap;
    elseif iquad==3
        ukmap(1:512,513:1024)=sukmap;
        tmmap(1:512,513:1024)=stmmap;
    else
        ukmap(513:1024,513:1024)=sukmap;
        tmmap(513:1024,513:1024)=stmmap;
    end
end
%%
%%% get masks %%%%
bigmask=maskdat.mask(ifield).bigmask;
nosrcmask=maskdat.mask(ifield).nosrc;
%% plot to see the map
figure
setwinsize(gcf,1500,500)
subplot(1,3,1)
imageclip(log10(abs(calmap.*cal)));
title('CIBER')
subplot(1,3,2)
imageclip(log10(abs(tmmap)));
title('2MASS sim map')
subplot(1,3,3)
imageclip(log10(abs(ukmap)));
title('UKIDSS sim map')
savename=strcat(savedir,'maps');
%print(savename,'-dpng');
%% Masking test: catalogs
alpha=-15;
m_max=17;
beta=-m_max*alpha+3;

pscmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','PSC');
xscmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','XSC');
xscrmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','XSCrej');
pscrmask1=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','A');
pscrmask2=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','B');
pscrmask3=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','C');
pscrmask4=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','D');
pscrmask5=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','E');
pscrmask6=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','F');
%%
figure
setwinsize(gcf,1000,1000)
subplot(3,3,1)
imageclip(pscmask);
title('PSC');
subplot(3,3,2)
imageclip(xscmask);
title('XSC');
subplot(3,3,3)
imageclip(xscrmask);
title('XSC rej');
subplot(3,3,4)
imageclip(pscrmask1);
title('PSC rej A');
subplot(3,3,5)
imageclip(pscrmask2);
title('PSC rej B');
subplot(3,3,6)
imageclip(pscrmask3);
title('PSC rej C');
subplot(3,3,7)
imageclip(pscrmask4);
title('PSC rej D');
subplot(3,3,8)
imageclip(pscrmask5);
title('PSC rej E');
subplot(3,3,9)
imageclip(pscrmask6);
title('PSC rej F');
savename=strcat(savedir,'masks');
print(savename,'-dpng');
%%
a=log10(ukmap);
b=a.*pscmask;
c=b.*xscmask.*xscrmask;
d=c.*pscrmask1;
e=d.*pscrmask2.*pscrmask3.*pscrmask4.*pscrmask5.*pscrmask6;
%% XSC demo
figure
setwinsize(gcf,1000,1000)
subplot(2,3,1)
imageclip(a);
title('unmasked');
v=caxis;
subplot(2,3,2)
imageclip(b);
title('PSC mask');
caxis(v);
subplot(2,3,3)
imageclip(c);
title('PSC & XSC mask');
caxis(v);
subplot(2,3,4)
imageclip(a(430:470,320:360));
caxis(v);
subplot(2,3,5)
imageclip(b(430:470,320:360));
caxis(v);
subplot(2,3,6)
imageclip(c(430:470,320:360));
caxis(v);
savename=strcat(savedir,'xsc_demo');
print(savename,'-dpng');
%% PSCrej demo
figure
setwinsize(gcf,1000,1000)
subplot(2,3,1)
imageclip(a);
title('unmasked');
v=caxis;
subplot(2,3,2)
imageclip(c);
title('PSC & XSC mask');
caxis(v);
subplot(2,3,3)
imageclip(d);
title('PSC & XSC & PSC rej A mask');
caxis(v);
subplot(2,3,4)
imageclip(a(400:430,580:610));
caxis(v);
subplot(2,3,5)
imageclip(c(400:430,580:610));
caxis(v);
subplot(2,3,6)
imageclip(d(400:430,580:610));
caxis(v);
savename=strcat(savedir,'pscr_demo');
print(savename,'-dpng');

%%
[~,~,~,~,binl]=get_angular_spec(randn(1024),randn(1024),pixscale);

figure
mask=nosrcmask.*pscmask;
[~,mask,~,~,~,clipmax,clipmin]=get_skymap(calmap.*cal,mask,5,5);
mask(find(ukmap>clipmax))=0;
ukmap1=(ukmap-mean(ukmap(find(mask)))).*mask;
[Cl,l]=get_angular_spec(ukmap1,ukmap1,pixscale);
maskfrac=numel(find(mask))./numel(mask);
mkk=numel(mask)./numel(find(mask));
Cl=Cl.*mkk;
%mkk=get_mkk_sim(mask,pixscale,binl,10,numel(binl),1,ones(1024),0,NaN);
%[Cl,wCl,~,l]=get_Cl(ukmap,mask,mkk,pixscale,ones(1024));
loglog(l,l.*(l+1).*Cl./2./pi,'-','Displayname','PSC mask');hold on

mask=nosrcmask.*pscmask.*xscmask.*xscrmask;
[~,mask,~,~,~,clipmax,clipmin]=get_skymap(calmap.*cal,mask,5,5);
mask(find(ukmap>clipmax))=0;
ukmap1=(ukmap-mean(ukmap(find(mask)))).*mask;
[Cl,l]=get_angular_spec(ukmap1,ukmap1,pixscale);
maskfrac=numel(find(mask))./numel(mask);
mkk=numel(mask)./numel(find(mask));
Cl=Cl.*mkk;
%mkk=get_mkk_sim(mask,pixscale,binl,10,numel(binl),1,ones(1024),0,NaN);
%[Cl,wCl,~,l]=get_Cl(ukmap,mask,mkk,pixscale,ones(1024));
loglog(l,l.*(l+1).*Cl./2./pi,'-','Displayname','PSC*XSC mask');hold on

mask=nosrcmask.*pscmask.*xscmask.*xscrmask;
mask=mask.*pscrmask1.*pscrmask2.*pscrmask3.*pscrmask4.*pscrmask5.*pscrmask6;
[~,mask,~,~,~,clipmax,clipmin]=get_skymap(calmap.*cal,mask,5,5);
mask(find(ukmap>clipmax))=0;
ukmap1=(ukmap-mean(ukmap(find(mask)))).*mask;
[Cl,l]=get_angular_spec(ukmap1,ukmap1,pixscale);
maskfrac=numel(find(mask))./numel(mask);
mkk=numel(mask)./numel(find(mask));
Cl=Cl.*mkk;
%mkk=get_mkk_sim(mask,pixscale,binl,10,numel(binl),1,ones(1024),0,NaN);
%[Cl,wCl,~,l]=get_Cl(ukmap,mask,mkk,pixscale,ones(1024));
loglog(l,l.*(l+1).*Cl./2./pi,'-','Displayname','PSC*XSC*all rej mask');hold on

xlim([1e2,2e5]);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
leg=legend('show','location','northwest');
legend boxoff
set(leg,'FontSize',15);
savename=strcat(savedir,'Clcats');
%print(savename,'-dpng');