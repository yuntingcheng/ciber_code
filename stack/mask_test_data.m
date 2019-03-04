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
%% get the bl from psfmap
loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');

load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');
ifield=4;

bestparam=fitpsfdat(ifield).bestparam_norm;
A=bestparam(1);
B=bestparam(2);
sig=bestparam(3);
r0=bestparam(4);
alpha=bestparam(5);
npix=512;

radmap = make_radius_map(zeros(2*npix+1),npix,npix).*pixsize;
psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);
[bl1,l]=get_angular_spec(psfmap,psfmap,pixscale);

radmap = make_radius_map(zeros(2*npix+1),npix+0.25,npix+0.25).*pixsize;
psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);
[bl2,l]=get_angular_spec(psfmap,psfmap,pixscale);

radmap = make_radius_map(zeros(2*npix+1),npix+0.5,npix+0.5).*pixsize;
psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);
[bl3,l]=get_angular_spec(psfmap,psfmap,pixscale);
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

%%% get sim srcmap %%%
ukmap=zeros(1024);
tmmap=zeros(1024);
for iquad=1:4
    quad=quad_arr(iquad);
    stmmap=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapt_2m.fits'));
    
    if iquad==1
        tmmap(1:512,1:512)=stmmap;
    elseif iquad==2
        tmmap(513:1024,1:512)=stmmap;
    elseif iquad==3
        tmmap(1:512,513:1024)=stmmap;
    else
        tmmap(513:1024,513:1024)=stmmap;
    end
end
%%% get masks %%%%
bigmask=maskdat.mask(ifield).bigmask;
nosrcmask=maskdat.mask(ifield).nosrc;
%%
alpha=-15;
m_arr=15.5:0.5:17;
for i=1:numel(m_arr)
    m_max=m_arr(i);
    pscmask1=make_strmask_2m(flight,inst,ifield,alpha,250,m_max,'catname','PSC');
    tmmap1=(tmmap-mean(tmmap(find(pscmask1)))).*pscmask1;
    [Cl1,l]=get_angular_spec(tmmap1,tmmap1,pixscale);
    mkk1=numel(pscmask1)./numel(find(pscmask1));
    %Cl1=Cl1.*mkk1;
    loglog(l,l.*(l+1).*Cl1./2./pi,'color',get_color(i));hold on
    drawnow
    pscmask1=make_strmask_2m(flight,inst,ifield,alpha,350,m_max,'catname','PSC');
    tmmap1=(tmmap-mean(tmmap(find(pscmask1)))).*pscmask1;
    [Cl1,l]=get_angular_spec(tmmap1,tmmap1,pixscale);
    mkk1=numel(pscmask1)./numel(find(pscmask1));
    %Cl1=Cl1.*mkk1;
    loglog(l,l.*(l+1).*Cl1./2./pi,'o--','color',get_color(i));hold on
    drawnow

end
%%
loglog(l,l.*(l+1).*bl1.*1e11./2./pi,'color','k','LineWidth',3);
loglog(l,l.*(l+1).*bl2.*1e11./2./pi,'color','k','LineWidth',3);
loglog(l,l.*(l+1).*bl3.*1e11./2./pi,'color','k','LineWidth',3);

%%
xlim([1e2,2e5]);
ylim([1e-2,1e4]);

xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
savename=strcat(savedir,'mask_test_dat');
print(savename,'-dpng');