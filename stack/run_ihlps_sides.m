flight=40030;
inst=1;
interp = 1;
mypaths=get_paths(flight);

savedir = strcat(mypaths.ciberdir, 'doc/20171018_stackihl/ihl_sides/plots/');
loaddir = strcat(mypaths.ciberdir, 'doc/20171018_stackihl/ihl_sides/TM',...
    num2str(inst),'/');
srcmapdir = strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
ifield=8;
dt=get_dark_times(flight,inst,ifield);
%% make the mask and MKK
pixscale=7;
[~,~,~,~,binl]=get_angular_spec(randn(720),randn(720),pixscale);
maskf = ones(720);
for m_min = 16:19
    m_max = m_min + 1;
    maskf = maskf .* fits_read(strcat(loaddir,dt.name,'_mask_frac',...
    num2str(m_min),'_',num2str(m_max),'.fits'));
end
mkkf =get_mkk_sim(maskf,pixscale,binl,10,numel(binl),1,ones(720),0,NaN);
mask = make_mask_sides(flight,inst,13,20);
mkk =get_mkk_sim(mask,pixscale,binl,10,numel(binl),1,ones(720),0,NaN);
%% get the sources map
% all the sources
psfmap_all = zeros(720);
ihlmap_all = zeros(720);
for m_min = 13:28
    m_max = m_min + 1;
    srcmap = fits_read(strcat(srcmapdir,dt.name,'_srcmap',...
    num2str(m_min),'_',num2str(m_max),'sides.fits'));
    ihlmap = fits_read(strcat(loaddir,dt.name,'_ihlmap_all',...
    num2str(m_min),'_',num2str(m_max),'.fits'));
    ihlmap(find(ihlmap==Inf))=0;
    psfmap_all = psfmap_all + srcmap;
    ihlmap_all = ihlmap_all  + ihlmap + srcmap;
end

% all the srouces m<20
psfmap_ps = zeros(720);
ihlmap_ps = zeros(720);
for m_min = 13:19
    m_max = m_min + 1;
    srcmap = fits_read(strcat(srcmapdir,dt.name,'_srcmap',...
    num2str(m_min),'_',num2str(m_max),'sides.fits'));
    ihlmap = fits_read(strcat(loaddir,dt.name,'_ihlmap_all',...
    num2str(m_min),'_',num2str(m_max),'.fits'));
    ihlmap(find(ihlmap==Inf))=0;
    psfmap_ps = psfmap_ps + srcmap;
    ihlmap_ps = ihlmap_ps + ihlmap + srcmap;
end

% frac of the srouces 16<m<20
psfmap_f = zeros(720);
ihlmap_f = zeros(720);
for m_min = 16:19
    m_max = m_min + 1;
    srcmap = fits_read(strcat(loaddir,dt.name,'_srcmap_frac',...
    num2str(m_min),'_',num2str(m_max),'.fits'));
    ihlmap = fits_read(strcat(loaddir,dt.name,'_ihlmap_frac',...
    num2str(m_min),'_',num2str(m_max),'.fits'));
    ihlmap(find(ihlmap==Inf))=0;
    psfmap_f = psfmap_f + srcmap;
    ihlmap_f = ihlmap_f + ihlmap + srcmap;
end
%% plot the sides PS
[Clihl_f,~,~,l]=get_Cl(ihlmap_f.*maskf,maskf,mkkf,pixscale,ones(720));
[Clpsf_f]=get_Cl(psfmap_f.*maskf,maskf,mkkf,pixscale,ones(720));

[Clihl_ps,~,~,l]=get_Cl(ihlmap_ps.*mask,mask,mkk,pixscale,ones(720));
[Clpsf_ps]=get_Cl(psfmap_ps.*mask,mask,mkk,pixscale,ones(720));

[Clihl_all,~,~,l]=get_Cl(ihlmap_all.*mask,mask,mkk,pixscale,ones(720));
[Clpsf_all]=get_Cl(psfmap_all.*mask,mask,mkk,pixscale,ones(720));

excessdir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/excessdat',excessdir),'excessdat');

excessdat.sides.l = l;
excessdat.sides.Clihl_f = Clihl_f;
excessdat.sides.Clpsf_f = Clpsf_f;
excessdat.sides.Clihl_ps = Clihl_ps;
excessdat.sides.Clpsf_ps = Clpsf_ps;
excessdat.sides.Clihl_all = Clihl_all;
excessdat.sides.Clpsf_all = Clpsf_all;

% get the corr factor for CIBER power spectrum
r_ps_sides = (Clihl_ps - Clpsf_ps)./ (Clihl_f - Clpsf_f);
r_all_sides = (Clihl_all - Clpsf_all)./ (Clihl_f - Clpsf_f);
l_cb = excessdat.Cldat(1).l;

r_ps = interp1(l,r_ps_sides,l_cb);
r_all = interp1(l,r_all_sides,l_cb);

excessdat.sides.r_ps = r_ps;
excessdat.sides.r_all = r_all;

save(sprintf('%s/excessdat',excessdir),'excessdat');

figure
setwinsize(gcf,800,300)

subplot(1,2,1)
loglog(l,l.*(l+1).*Clpsf_f./2./pi,'--','color',get_color(1),'linewidth',2,...
    'DisplayName','16 < m < 20 in catelog PSF');hold on
loglog(l,l.*(l+1).*Clihl_f./2./pi,'-','color',get_color(1),'linewidth',2,...
    'DisplayName','16 < m < 20 in catelog PSF + IHL');
loglog(l,l.*(l+1).*Clpsf_ps./2./pi,'--','color',get_color(2),'linewidth',2,...
    'DisplayName','m < 20 PSF');hold on
loglog(l,l.*(l+1).*Clihl_ps./2./pi,'-','color',get_color(2),'linewidth',2,...
    'DisplayName','m < 20 PSF + IHL');
loglog(l,l.*(l+1).*Clpsf_all./2./pi,'--','color',get_color(3),'linewidth',2,...
    'DisplayName','m < 29 PSF');hold on
loglog(l,l.*(l+1).*Clihl_all./2./pi,'-','color',get_color(3),'linewidth',2,...
    'DisplayName','m < 29 PSF + IHL');

xlim([1e2,2e5]);
ylim([1e-4,1e2]);
title(sprintf('Power Spectrum'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18);
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18);
h=legend('show','Location','northwest');
set(h,'fontsize',6)

subplot(1,2,2)
loglog(l,l.*(l+1).*(Clihl_f-Clpsf_f)./2./pi,'-','color',get_color(1),...
    'linewidth',2,'DisplayName','16 < m < 20 in catelog PSF + IHL');hold on
loglog(l,l.*(l+1).*(Clihl_ps-Clpsf_ps)./2./pi,'-','color',get_color(2),...
    'linewidth',2,'DisplayName','m < 20 Clall - Clpsf');
loglog(l,l.*(l+1).*(Clihl_all-Clpsf_all)./2./pi,'-','color',get_color(3),...
    'linewidth',2,'DisplayName','m < 29 Clall - Clpsf');

xlim([1e2,2e5]);
title(sprintf('Power Excess from IHL'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18);
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18);
h=legend('show','Location','southeast');
set(h,'fontsize',6)

savename=strcat(savedir,'IHL_power_spec_sides');
print(savename,'-dpng','-r250');%close
%%
flight=40030;
mypaths=get_paths(flight);

inst=1;
ifield=8;
dt=get_dark_times(flight,inst,ifield);

if inst == 1
    c_min = 0.96;
    c_mid = 1.05;
    c_max = 1.15;
elseif inst==2
    c_min = 1.59;
    c_mid = 2.45;
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
    
    subplot(2,2,1)
    loglog(l,l.*(l+1).*Clpsf0./2./pi,'--','color',get_color(im),'linewidth',2,...
        'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF'));hold on
    loglog(l,l.*(l+1).*Clihl0./2./pi,'-','color',get_color(im),'linewidth',2,...
        'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF+IHL'));
    
    subplot(2,2,2)
    loglog(l,l.*(l+1).*(Clihl0 - Clpsf0)./2./pi,'-','color',get_color(im),...
    'linewidth',2,'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),...
    ' Clall - Clpsf'));hold on
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
' Clall - Clpsf'));hold on

subplot(2,2,3)
r_ps = excessdat.sides.r_ps;
r_all = excessdat.sides.r_all;
loglog(l,l.*(l+1).*(Clihl0 - Clpsf0)./2./pi,'k--',...
'linewidth',2,'DisplayName',strcat('16<m<20 before correction'));hold on
loglog(l,l.*(l+1).*(Clihl0 - Clpsf0).*r_ps./2./pi,'r',...
   'linewidth',2,'DisplayName',...
strcat('16<m<20'));hold on
loglog(l,l.*(l+1).*(Clihl0 - Clpsf0).*r_all./2./pi,'b',...
    'linewidth',2,'DisplayName',...
 strcat('m<29'));hold on

subplot(2,2,4)
ymin = l.*(l+1).*(Clihl0 - Clpsf0).*c_min.^2.*r_ps./2./pi;
ymid = l.*(l+1).*(Clihl0 - Clpsf0).*c_mid.^2.*r_ps./2./pi;
ymax = l.*(l+1).*(Clihl0 - Clpsf0).*c_max.^2.*r_ps./2./pi;
y0 = l.*(l+1).*(Clihl0 - Clpsf0).*r_ps./2./pi;
loglog(l,y0,'b','linewidth',2,'DisplayName',...
 strcat('16<m<20'));hold on
loglog(l,ymin,'b:','linewidth',2,'DisplayName',...
 strcat('16<m<20'));hold on
loglog(l,ymax,'b:','linewidth',2,'DisplayName',...
 strcat('16<m<20'));hold on
ymin = l.*(l+1).*(Clihl0 - Clpsf0).*c_min.^2.*r_all./2./pi;
ymid = l.*(l+1).*(Clihl0 - Clpsf0).*c_mid.^2.*r_all./2./pi;
ymax = l.*(l+1).*(Clihl0 - Clpsf0).*c_max.^2.*r_all./2./pi;
y0 = l.*(l+1).*(Clihl0 - Clpsf0).*r_all./2./pi;
loglog(l,y0,'r','linewidth',2,'DisplayName',...
 strcat('m<20'));hold on
loglog(l,ymin,'r:','linewidth',2,'DisplayName',...
 strcat('m<20'));hold on
loglog(l,ymax,'r:','linewidth',2,'DisplayName',...
 strcat('m<20'));hold on

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
ylim([1e-3,1e0]);
title(sprintf('Color Interp Correction'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northeastoutside');
set(h,'fontsize',12)

savename=strcat(savedir,'IHL_power_spec');
print(savename,'-dpng','-r250');%close