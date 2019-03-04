flight=40030;
inst=2;
mypaths=get_paths(flight);
%%

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/excessdat',loaddir),'excessdat');

% fit two power law between 50-500 arcsec
fitmin = 8;
fitmid = 30;
fitmax = 200;

figure
setwinsize(gcf,1000,800)

for im=10:13
    m_min = excessdat.m_min_arr(im);
    m_max = excessdat.m_max_arr(im);
    r_arr = excessdat.r_arr;
    prof = excessdat.avg(im).prof;
    err = excessdat.avg(im).err;
    [fit,fit_model] = fit_pl2_prof_chi2(r_arr,prof,err,fitmin,fitmid,fitmax);

    excessdat.avg(im).fit_params = fit.params;
    excessdat.avg(im).fit_radius = fit.radius;
    excessdat.avg(im).fit_model = fit_model;
    
    subplot(2,2,im-9)
    loglog(r_arr,prof,'.','markersize',10);hold on
    errorbar(r_arr,prof,err,'k.');
    loglog(r_arr,-prof,'ko','markersize',5);
    errorbar(r_arr,-prof,err,'ko','markersize',5);
    h=loglog(r_arr,fit_model,'r','DisplayName',...
        sprintf('slope=(%.2f, %.2f)',fit.params(2),fit.params(3)));
    vline(fitmin,'b--');
    vline(fitmid,'b--');
    vline(fitmax,'b--');
    leg=legend(h,'Location','northeast');
    set(leg,'fontsize',10)
    legend boxoff

    xlim([4e-1,1e3])
    ylim([1e-4,1e0])
    xlabel('arcsec')
    ylabel('normalized <I_{stack}>')
    title(strcat(num2str(m_min),'<mAB(y band)<',num2str(m_max)));
end

pltsavedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));
savename=strcat(pltsavedir,'excess_pwl_fit');
print(savename,'-dpng');close

save(sprintf('%s/excessdat',loaddir),'excessdat');

%% make IHL srcmap

flight=40030;
inst=2;
interp = 1;
mypaths=get_paths(flight);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/excessdat',loaddir),'excessdat');

savedir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/ihl_pl/TM',...
    num2str(inst),'/');

ifield=8;
dt=get_dark_times(flight,inst,ifield);
for im=10:13
    m_min = excessdat.m_min_arr(im);
    m_max = excessdat.m_max_arr(im);
    params = excessdat.avg(im).fit_params;
    radius = excessdat.avg(im).fit_radius;
    ihlmap = make_ihl_pl2(flight,inst,ifield,m_min,m_max,params,radius,interp);
    fits_write(strcat(savedir,dt.name,'_ihlmap_pl',...
    num2str(m_min),'_',num2str(m_max)),ihlmap);
end
%%
ifield=8;
dt=get_dark_times(flight,inst,ifield);

srcmapdir = strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
savedir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/ihl_pl/TM',...
    num2str(inst),'/');

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
load(sprintf('%s/excessdat',loaddir),'excessdat');

countg_arr = stackmapdat(ifield).count_stackg_arr;

pixscale=7;
[~,~,~,~,binl]=get_angular_spec(randn(1024),randn(1024),pixscale);

figure
setwinsize(gcf,800,400)
Clpsftot = zeros(1,21);
Clihltot = zeros(1,21);

for im=10:13
    m_min = excessdat.m_min_arr(im);
    m_max = excessdat.m_max_arr(im);
        
    srcmap = fits_read(strcat(srcmapdir,dt.name,'_srcmapg',...
         num2str(m_min),'_',num2str(m_max),'ps.fits'));
    ihlmap = fits_read(strcat(savedir,dt.name,'_ihlmap_pl',...
    num2str(m_min),'_',num2str(m_max),'.fits'));
    ihlmap(find(ihlmap==Inf)) = 0;
    
    mask = make_mask_ps(flight,inst,ifield,1,m_min,m_max);
    mkk =get_mkk_sim(mask,pixscale,binl,10,numel(binl),1,ones(1024),0,NaN);
    
    mmap = (srcmap+ihlmap).*mask;
    [Clihl0,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
    
    mmap = (srcmap).*mask;
    [Clpsf0,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
    
    N_helg = IGLcounts_helgason(inst,(m_min+m_max)/2) * 4; % CIBER field 4 deg^2
    N_stack = countg_arr(im); % CIBER field 4 deg^2
    Nratio = N_helg/N_stack;
    Clihl = Clihl0.*Nratio.^2;
    Clpsf = Clpsf0.*Nratio.^2;
    
    Clpsftot = Clpsftot + Clpsf;
    Clihltot = Clihltot + Clihl;
    
    loglog(l,l.*(l+1).*Clpsf./2./pi,'--','color',get_color(im),'linewidth',2,...
        'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF'));hold on
    loglog(l,l.*(l+1).*Clihl./2./pi,'-','color',get_color(im),'linewidth',2,...
        'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF+IHL'));
    %loglog(l,l.*(l+1).*Clihl0./2./pi,':','color',get_color(im),'linewidth',2,...
    %    'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),...
    %    ' PSF+IHL (raw counts)'));
    drawnow
end
loglog(l,l.*(l+1).*Clpsftot./2./pi,'k--','linewidth',2,...
    'DisplayName','tot PSF');hold on
loglog(l,l.*(l+1).*Clihltot./2./pi,'k-','linewidth',2,...
    'DisplayName','tot PSF+IHL');

xlim([1e2,2e5]);
ylim([1e-4,1e0]);
title(sprintf('double power-law IHL'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northeastoutside');
set(h,'fontsize',12)

pltsavedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));
savename=strcat(pltsavedir,'IHL_power_spec');
print(savename,'-dpng');%close