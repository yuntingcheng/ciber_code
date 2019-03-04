flight=40030;
inst=1;
mypaths=get_paths(flight);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/excessdat',loaddir),'excessdat');

% fit two power law between 50-500 arcsec
fitmin = 8;
fitmid = 30;
fitmax = 200;

figure

m_min = excessdat.m_min_arr(10);
m_max = excessdat.m_max_arr(12);
r_arr = excessdat.r_arr;
prof = excessdat.all_avg.prof;
err = excessdat.all_avg.err;
[fit,fit_model] = fit_pl2_prof_chi2(r_arr,prof,err,fitmin,fitmid,fitmax);

excessdat.all_avg.fit_params = fit.params;
excessdat.all_avg.fit_radius = fit.radius;
excessdat.all_avg.fit_model = fit_model;

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


pltsavedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));
savename=strcat(pltsavedir,'excess_pwl_fit_all');
print(savename,'-dpng');%close

save(sprintf('%s/excessdat',loaddir),'excessdat');
%% make IHL srcmap
flight=40030;
interp = 1;
mypaths=get_paths(flight);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/excessdat',loaddir),'excessdat');

savedir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/ihl_pl/TM',...
    num2str(inst),'/');

ifield=8;
dt=get_dark_times(flight,inst,ifield);
for im=10:14
    m_min = excessdat.m_min_arr(im);
    m_max = excessdat.m_max_arr(im);
    params = excessdat.all_avg.fit_params;
    radius = excessdat.all_avg.fit_radius;
    ihlmap = make_ihl_pl2(flight,inst,ifield,m_min,m_max,params,radius,interp);
    fits_write(strcat(savedir,dt.name,'_ihlmap_pl',...
    num2str(m_min),'_',num2str(m_max)),ihlmap);
end
%% calculate and save PS
flight=40030;
mypaths=get_paths(flight);

inst=1;
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

psfmaptot = zeros(1024);
ihlmaptot = zeros(1024);
masktot = ones(1024);
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
    [Clihl,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
    
    mmap = (srcmap).*mask;
    [Clpsf,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
    
    N_helg = IGLcounts_helgason(inst,(m_min+m_max)/2) * 4; % CIBER field 4 deg^2
    N_stack = countg_arr(im); % CIBER field 4 deg^2
    Nratio = N_helg/N_stack;
    
    Cldat(im).l = l;
    Cldat(im).Clpsf = Clpsf;
    Cldat(im).Clihl = Clihl;
    Cldat(im).Nratio = Nratio;
    
    psfmaptot = psfmaptot + srcmap;
    ihlmaptot = ihlmaptot + ihlmap;
    masktot = masktot.*mask;
end

mkk = get_mkk_sim(masktot,pixscale,binl,10,numel(binl),1,ones(1024),0,NaN);
mmap = (psfmaptot + ihlmaptot).*masktot;
[Clihl,~,~,l]=get_Cl(mmap,masktot,mkk,pixscale,ones(1024));

mmap = psfmaptot.*masktot;
[Clpsf,~,~,l]=get_Cl(mmap,masktot,mkk,pixscale,ones(1024));

Cldat(1).l = l;
Cldat(1).Clpsf = Clpsf;
Cldat(1).Clihl = Clihl;

excessdat.Cldat = Cldat;
save(sprintf('%s/excessdat',loaddir),'excessdat');
%% plot results

flight=40030;
mypaths=get_paths(flight);

inst=2;
ifield=8;
dt=get_dark_times(flight,inst,ifield);

if inst == 1
    c_min = 0.88;
    c_mid = 1.34;
    c_max = 1.81;
elseif inst==2
    c_min = 1.20;
    c_mid = 2.25;
    c_max = 3.30;
end

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/excessdat',loaddir),'excessdat');

figure
setwinsize(gcf,1500,600)
Nratio_arr = [];

for im=10:13
    m_min = excessdat.m_min_arr(im);
    m_max = excessdat.m_max_arr(im);
    
    l = excessdat.Cldat(im).l;
    Clpsf0 = excessdat.Cldat(im).Clpsf;
    Clihl0 = excessdat.Cldat(im).Clihl;
    Nratio = excessdat.Cldat(im).Nratio;
    Nratio_arr = [Nratio Nratio_arr];
    
    subplot(2,2,1)
    loglog(l,l.*(l+1).*Clpsf0./2./pi,'--','color',get_color(im),'linewidth',2,...
        'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF'));hold on
    loglog(l,l.*(l+1).*Clihl0./2./pi,'-','color',get_color(im),'linewidth',2,...
        'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF+IHL'));
    
    subplot(2,2,2)
    loglog(l,l.*(l+1).*(Clihl0 - Clpsf0)./2./pi,'-','color',get_color(im),...
    'linewidth',2,'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),...
    'Clpsf - Clall'));hold on
    
    subplot(2,2,3)
    loglog(l,l.*(l+1).*(Clihl0 - Clpsf0).*Nratio./2./pi,'--',...
        'color',get_color(im),'linewidth',2,'DisplayName',...
   strcat(num2str(m_min),'<m<',num2str(m_max),'(Clpsf - Clall)Nr'));hold on
    loglog(l,l.*(l+1).*(Clihl0 - Clpsf0).*Nratio.^2./2./pi,'-',...
        'color',get_color(im),'linewidth',2,'DisplayName',...
     strcat(num2str(m_min),'<m<',num2str(m_max),'(Clpsf - Clall)Nr^2'));hold on    
    
    drawnow
end

m_min = excessdat.m_min_arr(10);
m_max = excessdat.m_max_arr(13);

l = excessdat.Cldat(1).l;
Clpsf0 = excessdat.Cldat(1).Clpsf;
Clihl0 = excessdat.Cldat(1).Clihl;

subplot(2,2,1)
loglog(l,l.*(l+1).*Clpsf0./2./pi,'k--','linewidth',2,...
    'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF'));hold on
loglog(l,l.*(l+1).*Clihl0./2./pi,'k-','linewidth',2,...
    'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF+IHL'));
subplot(2,2,2)
loglog(l,l.*(l+1).*(Clihl0 - Clpsf0)./2./pi,'k-',...
'linewidth',2,'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),...
'Clpsf - Clall'));hold on

subplot(2,2,3)
loglog(l,l.*(l+1).*(Clihl0 - Clpsf0).*min(Nratio_arr)./2./pi,'k--',...
   'linewidth',2,'DisplayName',...
strcat(num2str(m_min),'<m<',num2str(m_max),'(Clpsf - Clall)Nr'));hold on
loglog(l,l.*(l+1).*(Clihl0 - Clpsf0).*max(Nratio_arr).^2./2./pi,'k-',...
    'linewidth',2,'DisplayName',...
 strcat(num2str(m_min),'<m<',num2str(m_max),'(Clpsf - Clall)Nr^2'));hold on

subplot(2,2,4)
ymin = l.*(l+1).*(Clihl0 - Clpsf0).*c_min.*2.*max(Nratio_arr).^2./2./pi;
ymid = l.*(l+1).*(Clihl0 - Clpsf0).*c_mid.*2.*max(Nratio_arr).^2./2./pi;
ymax = l.*(l+1).*(Clihl0 - Clpsf0).*c_max.*2.*max(Nratio_arr).^2./2./pi;
y0 = l.*(l+1).*(Clihl0 - Clpsf0).*max(Nratio_arr).^2./2./pi;
loglog(l,y0,'k--','linewidth',2,'DisplayName',...
 strcat(num2str(m_min),'<m<',num2str(m_max),'before correction'));hold on
%loglog(l,ymid,'k-','linewidth',2,'DisplayName',...
% strcat(num2str(m_min),'<m<',num2str(m_max),'after correction'));hold on
loglog(l,ymin,'k:','linewidth',2,'DisplayName',...
 strcat(num2str(m_min),'<m<',num2str(m_max),'after correction'));hold on
loglog(l,ymax,'k:','linewidth',2,'DisplayName',...
 strcat(num2str(m_min),'<m<',num2str(m_max),'after correction'));hold on

subplot(2,2,1)
xlim([1e2,2e5]);
ylim([1e-4,1e0]);
title(sprintf('Power Spectrum'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northeastoutside');
set(h,'fontsize',12)

subplot(2,2,2)
xlim([1e2,2e5]);
ylim([1e-4,1e0]);
title(sprintf('Power Excess from IHL'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northeastoutside');
set(h,'fontsize',12)

subplot(2,2,3)
xlim([1e2,2e5]);
ylim([1e-4,1e0]);
title(sprintf('Number Counts Correction'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northeastoutside');
set(h,'fontsize',12)

subplot(2,2,4)
xlim([1e2,2e5]);
ylim([1e-2,1e0]);
title(sprintf('Color Interp Correction (Nr^2)'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northeastoutside');
set(h,'fontsize',12)

pltsavedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));
savename=strcat(pltsavedir,'IHL_power_spec_all');
print(savename,'-dpng');%close
