flight=40030;
inst=1;
mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

fitmin = 7;
fitmax = inf;

figure
setwinsize(gcf,1200,600)
for im = 10:12

prof = 0;
err = 0;
for ifield = 4:8
    dt=get_dark_times(flight,inst,ifield);
    loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
    load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');
    diff = ihlprofdat.excess(im).diff;
    diff_err = ihlprofdat.excess(im).diff_err;
    prof = prof + diff./diff_err.^2;
    err = err + 1./diff_err.^2;
end
prof = prof ./ err;
err = sqrt(1./err);

r_arr = ihlprofdat.r_arr;
m_min = ihlprofdat.m_min_arr(im);
m_max = ihlprofdat.m_max_arr(im);


[fit] = fit_beta_model(r_arr,prof,err,fitmin,fitmax);

ihlsimdat(im).m_min = m_min;
ihlsimdat(im).m_max = m_max;
ihlsimdat(im).param = fit.params;

subplot(2,3,im-9)
loglog(r_arr,prof,'.','markersize',10);hold on
errorbar(r_arr,prof,err,'k.');
loglog(r_arr,-prof,'ko','markersize',5);
errorbar(r_arr,-prof,err,'ko','markersize',5);
h(1)=loglog(r_arr,fit.fit_model,'r','linewidth',2);
h(2) = vline(fitmin,'b');
vline(fitmax,'b');
leg=legend(h,{sprintf('rc = %.2f\nbeta = %.2f',fit.params(2),fit.params(3)),...
    'fitting range'},'Location','northeast');
set(leg,'fontsize',15)
legend boxoff
title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15)
xlim([4e-1,1e3])
ylim([1e-4,1e0])
xlabel('arcsec', 'fontsize',20);
ylabel('normalized <I_{stack}>','fontsize',20);

subplot(2,3,im-6)
semilogx(r_arr,prof,'.','markersize',10);hold on
errorbar(r_arr,prof,err,'k.');
Ylim = get(gca,'YLim');
h(1)=semilogx(r_arr,fit.fit_model,'r','linewidth',2);
h(2) = vline(fitmin,'b');
vline(fitmax,'b');
title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15)
xlim([4e-1,1e3])
ylim(Ylim)
xlabel('arcsec', 'fontsize',20);
ylabel('normalized <I_{stack}>','fontsize',20);

end
savename=strcat(pltsavedir,'rprof_excess_avg_fit');
print(savename,'-dpng');%close
%% run ihl map
savedir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
ifield = 8;
interp = 1;
for im = 10:12
m_min = ihlprofdat.m_min_arr(im);
m_max = ihlprofdat.m_max_arr(im);
params = ihlsimdat(im).param;

ihlmap = make_ihl_beta(flight,inst,ifield,m_min,m_max,params,inf,interp);
fits_write(strcat(savedir,dt.name,'_ihlmap_pl',...
num2str(m_min),'_',num2str(m_max)),ihlmap);

ihlmap = make_ihl_beta(flight,inst,ifield,m_min,m_max,params,100,interp);
fits_write(strcat(savedir,dt.name,'_ihlmap_pl',...
num2str(m_min),'_',num2str(m_max),'_100'),ihlmap);

end
%%
srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
ihlmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');

pixscale=7;
[~,~,~,~,binl]=get_angular_spec(randn(1024),randn(1024),pixscale);

psfmaptot = zeros(1024);
ihlmaptot = zeros(1024);
ihlmaptot100 = zeros(1024);
masktot = ones(1024);
for im = 10:12
    m_min = ihlsimdat(im).m_min;
    m_max = ihlsimdat(im).m_max;
    psfmap = fits_read(strcat(srcmapdir,dt.name,'_srcmapg',...
        num2str(m_min),'_',num2str(m_max),'ps.fits'));
    ihlmap = fits_read(strcat(ihlmapdir,dt.name,'_ihlmap_pl',...
        num2str(m_min),'_',num2str(m_max),'.fits'));
    ihlmap100 = fits_read(strcat(ihlmapdir,dt.name,'_ihlmap_pl',...
        num2str(m_min),'_',num2str(m_max),'_100.fits'));
    mask = make_mask_ps(flight,inst,ifield,1,m_min,m_max);
    mkk =get_mkk_sim(mask,pixscale,binl,10,numel(binl),1,ones(1024),0,NaN);
    ihlsimdat(im).psfmap = psfmap;
    ihlsimdat(im).ihlmap = ihlmap;
    ihlsimdat(im).ihlmap100 = ihlmap100;
    ihlsimdat(im).mask = mask;
    ihlsimdat(im).mkk = mkk;

    psfmaptot = psfmaptot + psfmap;
    ihlmaptot = ihlmaptot + ihlmap;
    ihlmaptot100 = ihlmaptot100 + ihlmap100;
    masktot = masktot.*mask;
end

mkk =get_mkk_sim(masktot,pixscale,binl,10,numel(binl),1,ones(1024),0,NaN);
ihlsimdat(1).psfmap = psfmaptot;
ihlsimdat(1).ihlmap = ihlmaptot;
ihlsimdat(1).ihlmap100 = ihlmaptot100;
ihlsimdat(1).mask = masktot;
ihlsimdat(1).mkk = mkk;
%%

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
countg_arr = stackmapdat(ifield).count_stackg_arr;

Nratio_arr = [];
figure
setwinsize(gcf,1500,400)
for im = 10:12
    m_min = ihlsimdat(im).m_min;
    m_max = ihlsimdat(im).m_max;
    psfmap = ihlsimdat(im).psfmap;
    ihlmap = ihlsimdat(im).ihlmap;
    ihlmap100 = ihlsimdat(im).ihlmap100;
    mask = ihlsimdat(im).mask;
    mkk = ihlsimdat(im).mkk;
    
    mmap = (psfmap).*mask;
    [Clpsf,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
    mmap = (psfmap+ihlmap).*mask;
    [Clihl,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
    mmap = (psfmap+ihlmap100).*mask;
    [Clihl100,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
    
    N_helg = IGLcounts_helgason(inst,(m_min+m_max)/2) * 4; % CIBER field 4 deg^2
    N_stack = countg_arr(im); % CIBER field 4 deg^2
    Nratio = N_helg/N_stack;
    Nratio_arr = [Nratio_arr Nratio];
    subplot(1,3,1)
    loglog(l,l.*(l+1).*Clpsf./2./pi,'--','color',get_color(im-9),'linewidth',2,...
        'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF'));hold on
    loglog(l,l.*(l+1).*Clihl100./2./pi,':','color',get_color(im-9),'linewidth',2,...
        'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF+IHL(<100)'));
    loglog(l,l.*(l+1).*Clihl./2./pi,'-','color',get_color(im-9),'linewidth',2,...
        'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF+IHL'));
    subplot(1,3,2)
    loglog(l,l.*(l+1).*(Clihl-Clpsf)./2./pi,'-','color',get_color(im-9),...
        'linewidth',2,'DisplayName',...
        strcat(num2str(m_min),'<m<',num2str(m_max),' (PSF+IHL) - PSF'));hold on
    subplot(1,3,3)
    loglog(l,l.*(l+1).*(Clihl-Clpsf).*(Nratio^2)./2./pi,'-',...
        'color',get_color(im-9),'linewidth',2,'DisplayName',...
        strcat(num2str(m_min),'<m<',num2str(m_max),'N^2 correction'));hold on
    loglog(l,l.*(l+1).*(Clihl-Clpsf).*(Nratio)./2./pi,'--',...
        'color',get_color(im-9),'linewidth',2,'DisplayName',...
        strcat(num2str(m_min),'<m<',num2str(m_max),'N correction'));hold on

end

mkk = ihlsimdat(1).mkk;
psfmap = ihlsimdat(1).psfmap;
ihlmap = ihlsimdat(1).ihlmap;
ihlmap100 = ihlsimdat(1).ihlmap100;
mask = ihlsimdat(1).mask;
mmap = (psfmap).*mask;
[Clpsf,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
mmap = (psfmap+ihlmap).*mask;
[Clihl,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
mmap = (psfmap+ihlmap100).*mask;
[Clihl100,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));

subplot(1,3,1)
loglog(l,l.*(l+1).*Clpsf./2./pi,'k--','linewidth',2,...
    'DisplayName',strcat('16<m<19',' PSF'));hold on
loglog(l,l.*(l+1).*Clihl100./2./pi,'k:','linewidth',2,...
    'DisplayName',strcat('16<m<19',' PSF+IHL'));
loglog(l,l.*(l+1).*Clihl./2./pi,'k-','linewidth',2,...
    'DisplayName',strcat('16<m<19',' PSF+IHL'));
xlim([1e2,2e5]);
ylim([1e-5,1e-1]);
title(sprintf('Power Spectrum'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northwest');
set(h,'fontsize',7)

subplot(1,3,2)
loglog(l,l.*(l+1).*(Clihl-Clpsf)./2./pi,'k-',...
    'linewidth',2,'DisplayName',strcat('16<m<19',' (PSF+IHL) - PSF'));
xlim([1e2,2e5]);
ylim([1e-5,1e-1]);
title(sprintf('Power Excess from IHL'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northwest');
set(h,'fontsize',7)

subplot(1,3,3)
loglog(l,l.*(l+1).*(Clihl-Clpsf).*(max(Nratio_arr)^2)./2./pi,'k-',...
    'linewidth',2,'DisplayName',...
    strcat('16<m<19',' N^2 correction'));hold on
loglog(l,l.*(l+1).*(Clihl-Clpsf).*(min(Nratio_arr)^2)./2./pi,'k--',...
    'linewidth',2,'DisplayName',...
    strcat('16<m<19',' N correction'));hold on

xlim([1e2,2e5]);
ylim([1e-5,1e-1]);
title(sprintf('Number Counts Correction'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','southeast');
set(h,'fontsize',7)
savename=strcat(pltsavedir,'IHL_power_spec');
print(savename,'-dpng');%close
