%% make sim IHL, src map for completeness correction
%{
flight=40030;
inst=1;
mypaths=get_paths(flight);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/excessdat',loaddir),'excessdat');

savedir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/corr_sim/TM',...
    num2str(inst),'/maps/');

ifield=8;
dt=get_dark_times(flight,inst,ifield);
for im=10:13
    m_min = excessdat.m_min_arr(im);
    m_max = excessdat.m_max_arr(im);

    [map_all1,map_use1] = make_srcmap_sim(flight,inst,ifield,1,m_min,m_max);
    [map_all0,map_use0] = make_srcmap_sim(flight,inst,ifield,0,m_min,m_max);

    
    fits_write(strcat(savedir,dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'corr_all'),map_all1);
    fits_write(strcat(savedir,dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'corr_use'),map_use1);
    fits_write(strcat(savedir,dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'rand_all'),map_all0);
    fits_write(strcat(savedir,dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'rand_use'),map_use0);

    params = excessdat.avg(im).fit_params;
    radius = excessdat.avg(im).fit_radius;

    [map_all1,map_use1] = make_ihl_pl2_sim(flight,inst,ifield,1,...
    m_min,m_max,params,radius);
    [map_all0,map_use0] = make_ihl_pl2_sim(flight,inst,ifield,0,...
    m_min,m_max,params,radius);

    fits_write(strcat(savedir,dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'corr_all'),map_all1);
    fits_write(strcat(savedir,dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'corr_use'),map_use1);
    fits_write(strcat(savedir,dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'rand_all'),map_all0);
    fits_write(strcat(savedir,dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'rand_use'),map_use0);

end
%}
%% get the mask & mkk
%{
flight=40030;
inst=1;
mypaths=get_paths(flight);
savedir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/corr_sim/TM',...
    num2str(inst),'/maps/');

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/excessdat',loaddir),'excessdat');

ifield=8;
pixscale=7;
[~,~,~,~,binl]=get_angular_spec(randn(1024),randn(1024),pixscale);
dt=get_dark_times(flight,inst,ifield);


maskdat(1).corr.mask_all = ones(1024);
maskdat(1).corr.mask_use = ones(1024);
maskdat(1).rand.mask_all = ones(1024);
maskdat(1).rand.mask_use = ones(1024);

for im=10:13
    m_min = excessdat.m_min_arr(im);
    m_max = excessdat.m_max_arr(im);

    disp(sprintf('=========corr mask %d<m<%d=============',m_min,m_max));

    [mask_all,mask_use] = make_ihl_sim_mask...
    (flight,inst,ifield,1,m_min,m_max);
    mkk_all=get_mkk_sim(mask_all,7,binl,10,numel(binl),1,ones(1024),0,NaN);
    mkk_use=get_mkk_sim(mask_use,7,binl,10,numel(binl),1,ones(1024),0,NaN);
    maskdat(im).corr.mask_all = mask_all;
    maskdat(im).corr.mask_use = mask_use;
    maskdat(im).corr.mkk_all = mkk_all;
    maskdat(im).corr.mkk_use = mkk_use;
    maskdat(1).corr.mask_all = maskdat(1).corr.mask_all.*mask_all;
    maskdat(1).corr.mask_use = maskdat(1).corr.mask_use.*mask_use;

    disp(sprintf('=========rand mask %d<m<%d=============',m_min,m_max));
    [mask_all,mask_use] = make_ihl_sim_mask...
    (flight,inst,ifield,0,m_min,m_max);
    mkk_all=get_mkk_sim(mask_all,7,binl,10,numel(binl),1,ones(1024),0,NaN);
    mkk_use=get_mkk_sim(mask_use,7,binl,10,numel(binl),1,ones(1024),0,NaN);
    maskdat(im).rand.mask_all = mask_all;
    maskdat(im).rand.mask_use = mask_use;
    maskdat(im).rand.mkk_all = mkk_all;
    maskdat(im).rand.mkk_use = mkk_use;
    maskdat(1).rand.mask_all = maskdat(1).rand.mask_all.*mask_all;
    maskdat(1).rand.mask_use = maskdat(1).rand.mask_use.*mask_use;

end
mask_all = maskdat(1).corr.mask_all;
mask_use = maskdat(1).corr.mask_use;
mkk_all=get_mkk_sim(mask_all,7,binl,10,numel(binl),1,ones(1024),0,NaN);
mkk_use=get_mkk_sim(mask_use,7,binl,10,numel(binl),1,ones(1024),0,NaN);
maskdat(1).corr.mkk_all = mkk_all;
maskdat(1).corr.mkk_use = mkk_use;

mask_all = maskdat(1).rand.mask_all;
mask_use = maskdat(1).rand.mask_use;
mkk_all=get_mkk_sim(mask_all,7,binl,10,numel(binl),1,ones(1024),0,NaN);
mkk_use=get_mkk_sim(mask_use,7,binl,10,numel(binl),1,ones(1024),0,NaN);
maskdat(1).rand.mkk_all = mkk_all;
maskdat(1).rand.mkk_use = mkk_use;

save(sprintf('%s/maskdat',savedir),'maskdat');
%}
%% plot the maps
%{
flight=40030;
inst=1;
mypaths=get_paths(flight);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

savedir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/corr_sim/TM',...
    num2str(inst),'/');
load(sprintf('%s/maps/maskdat',savedir),'maskdat');

ifield=8;
dt=get_dark_times(flight,inst,ifield);
countg_arr = stackmapdat(ifield).count_stackg_arr;
totcorr_all = zeros(1024);
totcorr_use = zeros(1024);
totrand_all = zeros(1024);
totrand_use = zeros(1024);

for im=10:13
    figure
    setwinsize(gcf,1000,800)

    m_min = stackmapdat(ifield).m_min_arr(im);
    m_max = stackmapdat(ifield).m_max_arr(im);

    N_helg = IGLcounts_helgason(inst,(m_min+m_max)/2) * 4;
    N_stack = countg_arr(im);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    src_all = fits_read(strcat(savedir,'maps/',dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'rand_all.fits'));
    src_use = fits_read(strcat(savedir,'maps/',dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'rand_use.fits'));

    ihl_all = fits_read(strcat(savedir,'maps/',dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'rand_all.fits'));
    ihl_use = fits_read(strcat(savedir,'maps/',dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'rand_use.fits'));
    ihl_all(find(ihl_all==inf))=0;
    ihl_use(find(ihl_use==inf))=0;
    
    totrand_all = totrand_all + src_all + ihl_all;
    totrand_use = totrand_use + src_use + ihl_use;
    
    subplot(2,2,1)
    imageclip(log10(src_all + ihl_all));
    v=caxis;
    title(sprintf('%d<m<%d, random',m_min,m_max));
    subplot(2,2,3)
    imageclip(log10(src_use + ihl_use));
    caxis(v);
    title(sprintf('%.2f %% sources in stacking',N_stack/N_helg*100));

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    src_all = fits_read(strcat(savedir,'maps/',dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'corr_all.fits'));
    src_use = fits_read(strcat(savedir,'maps/',dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'corr_use.fits'));

    ihl_all = fits_read(strcat(savedir,'maps/',dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'corr_all.fits'));
    ihl_use = fits_read(strcat(savedir,'maps/',dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'corr_use.fits'));
    ihl_all(find(ihl_all==inf))=0;
    ihl_use(find(ihl_use==inf))=0;
    
    totcorr_all = totcorr_all + src_all + ihl_all;
    totcorr_use = totcorr_use + src_use + ihl_use;
    
    subplot(2,2,2)
    imageclip(log10(src_all + ihl_all));
    caxis(v);
    title(sprintf('%d<m<%d, correlated',m_min,m_max));
    subplot(2,2,4)
    imageclip(log10(src_use + ihl_use));
    caxis(v);
    title(sprintf('%.2f %% sources in stacking',N_stack/N_helg*100));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    suptitle(sprintf('%d<mAB(y band)<%d (%d total, %d in stacking)',...
        m_min,m_max,round(N_helg),N_stack));

    savename=strcat(savedir,dt.name,'_simmap',num2str(m_min),'_',num2str(m_max));
    print(savename,'-dpng');%close
end

figure
setwinsize(gcf,1000,800)

subplot(2,2,1)
imageclip(log10(totrand_all));
v=caxis;
title(sprintf('total, random'));
subplot(2,2,3)
imageclip(log10(totrand_use));
caxis(v);
title(sprintf('sources in stacking'));

subplot(2,2,2)
imageclip(log10(totcorr_all));
caxis(v);
title(sprintf('total, correlated'));
subplot(2,2,4)
imageclip(log10(totcorr_use));
caxis(v);
title(sprintf('sources in stacking'));

suptitle(strcat(num2str(16),'<mAB(y band)<',num2str(20)));

savename=strcat(savedir,dt.name,'_simmap_all');
print(savename,'-dpng');%close
%}
%% plot the power spec
%{
flight=40030;
inst=1;
mypaths=get_paths(flight);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

savedir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/corr_sim/TM',...
    num2str(inst),'/');
load(sprintf('%s/maps/maskdat',savedir),'maskdat');

ifield=8;
dt=get_dark_times(flight,inst,ifield);
countg_arr = stackmapdat(ifield).count_stackg_arr;

tot_all = zeros(1024);
tot_use = zeros(1024);

figure
setwinsize(gcf,800,400)

for im=10:13

    m_min = stackmapdat(ifield).m_min_arr(im);
    m_max = stackmapdat(ifield).m_max_arr(im);

    N_helg = IGLcounts_helgason(inst,(m_min+m_max)/2) * 4;
    N_stack = countg_arr(im);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    src_all = fits_read(strcat(savedir,'maps/',dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'rand_all.fits'));
    src_use = fits_read(strcat(savedir,'maps/',dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'rand_use.fits'));
    ihl_all = fits_read(strcat(savedir,'maps/',dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'rand_all.fits'));
    ihl_use = fits_read(strcat(savedir,'maps/',dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'rand_use.fits'));
    ihl_all(find(ihl_all==inf))=0;
    ihl_use(find(ihl_use==inf))=0;
    mask_all = maskdat(im).rand.mask_all;
    mask_use = maskdat(im).rand.mask_use;
    mkk_all = maskdat(im).rand.mkk_all;
    mkk_use = maskdat(im).rand.mkk_use;
    tot_all = tot_all + src_all + ihl_all;
    tot_use = tot_use + src_use + ihl_use;
    
    [Clall,~,~,l]=get_Cl(src_all+ihl_all,mask_all,mkk_all,7,ones(1024));
    [Cluse,~,~,l]=get_Cl(src_use+ihl_use,mask_use,mkk_use,7,ones(1024));
    
    ratiodat(im).l = l;
    ratiodat(im).rand = Clall./Cluse;
    
   loglog(l,l.*(l+1).*Clall./2./pi,'-','color',get_color(im-9),'linewidth',2,...
     'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' all'));hold on
   loglog(l,l.*(l+1).*Cluse./2./pi,'--','color',get_color(im-9),'linewidth',2,...
        'DisplayName',sprintf('%d/%d (%.2f %%) in stack',...
        N_stack,round(N_helg),N_stack/N_helg*100));
 
end
mask_all = maskdat(1).rand.mask_all;
mask_use = maskdat(1).rand.mask_use;
mkk_all = maskdat(1).rand.mkk_all;
mkk_use = maskdat(1).rand.mkk_use;

[Clall,~,~,l]=get_Cl(tot_all,mask_all,mkk_all,7,ones(1024));
[Cluse,~,~,l]=get_Cl(tot_use,mask_use,mkk_use,7,ones(1024));
ratiodat(1).l = l;
ratiodat(1).rand = Clall./Cluse;
loglog(l,l.*(l+1).*Clall./2./pi,'k-','linewidth',2,...
 'DisplayName','all');hold on
loglog(l,l.*(l+1).*Cluse./2./pi,'k--','linewidth',2,...
    'DisplayName',sprintf('in stack'));

xlim([1e2,2e5]);
ylim([1e-6,1e0]);
title(sprintf('random sources(double power-law IHL)'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northeastoutside');
set(h,'fontsize',12)
savename=strcat(savedir,dt.name,'_PS_rand');
print(savename,'-dpng');%close

%%%%%%%%%%%%%%%%%%%%%
tot_all = zeros(1024);
tot_use = zeros(1024);

figure
setwinsize(gcf,800,400)

for im=10:13

    m_min = stackmapdat(ifield).m_min_arr(im);
    m_max = stackmapdat(ifield).m_max_arr(im);

    N_helg = IGLcounts_helgason(inst,(m_min+m_max)/2) * 4;
    N_stack = countg_arr(im);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    src_all = fits_read(strcat(savedir,'maps/',dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'corr_all.fits'));
    src_use = fits_read(strcat(savedir,'maps/',dt.name,'_src',...
    num2str(m_min),'_',num2str(m_max),'corr_use.fits'));
    ihl_all = fits_read(strcat(savedir,'maps/',dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'corr_all.fits'));
    ihl_use = fits_read(strcat(savedir,'maps/',dt.name,'_ihl',...
    num2str(m_min),'_',num2str(m_max),'corr_use.fits'));
    ihl_all(find(ihl_all==inf))=0;
    ihl_use(find(ihl_use==inf))=0;
    mask_all = maskdat(im).corr.mask_all;
    mask_use = maskdat(im).corr.mask_use;
    mkk_all = maskdat(im).corr.mkk_all;
    mkk_use = maskdat(im).corr.mkk_use;
    tot_all = tot_all + src_all + ihl_all;
    tot_use = tot_use + src_use + ihl_use;
    
    [Clall,~,~,l]=get_Cl(src_all+ihl_all,mask_all,mkk_all,7,ones(1024));
    [Cluse,~,~,l]=get_Cl(src_use+ihl_use,mask_use,mkk_use,7,ones(1024));
   
    ratiodat(im).corr = Clall./Cluse;

   loglog(l,l.*(l+1).*Clall./2./pi,'-','color',get_color(im-9),'linewidth',2,...
     'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' all'));hold on
   loglog(l,l.*(l+1).*Cluse./2./pi,'--','color',get_color(im-9),'linewidth',2,...
        'DisplayName',sprintf('%d/%d (%.2f %%) in stack',...
        N_stack,round(N_helg),N_stack/N_helg*100));
 
end
mask_all = maskdat(1).corr.mask_all;
mask_use = maskdat(1).corr.mask_use;
mkk_all = maskdat(1).corr.mkk_all;
mkk_use = maskdat(1).corr.mkk_use;

[Clall,~,~,l]=get_Cl(tot_all,mask_all,mkk_all,7,ones(1024));
[Cluse,~,~,l]=get_Cl(tot_use,mask_use,mkk_use,7,ones(1024));
ratiodat(1).corr = Clall./Cluse;
loglog(l,l.*(l+1).*Clall./2./pi,'k-','linewidth',2,...
 'DisplayName','all');hold on
loglog(l,l.*(l+1).*Cluse./2./pi,'k--','linewidth',2,...
    'DisplayName',sprintf('in stack'));

xlim([1e2,2e5]);
ylim([1e-6,1e0]);
title(sprintf('correlated sources(double power-law IHL)'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northeastoutside');
set(h,'fontsize',12)
savename=strcat(savedir,dt.name,'_PS_corr');
print(savename,'-dpng');%close

save(sprintf('%s/maps/ratiodat',savedir),'ratiodat');
%}
%% plot the correction factor
%{
flight=40030;
inst=1;
ifield=8;
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
countg_arr = stackmapdat(ifield).count_stackg_arr;

savedir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/corr_sim/TM',...
    num2str(inst),'/');
load(sprintf('%s/maps/ratiodat',savedir),'ratiodat');

l = ratiodat(1).l;
figure
setwinsize(gcf,800,600)

for im=10:13

subplot(2,2,im-9)
m_min = stackmapdat(ifield).m_min_arr(im);
m_max = stackmapdat(ifield).m_max_arr(im);
N_helg = IGLcounts_helgason(inst,(m_min+m_max)/2) * 4;
N_stack = countg_arr(im);
h(1) = semilogx(l, ratiodat(im).corr,'r-','linewidth',3);hold on
h(2) = semilogx(l, ratiodat(im).rand,'b-','linewidth',3);hold on
h(3) = hline((N_helg/N_stack),'k--');
h(4) = hline((N_helg/N_stack).^2,'k:');
h(3).LineWidth = 2;
h(4).LineWidth = 2;
xlim([1e2,2e5]);
title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$C_\ell(\rm all)/C_\ell(stacked)$',...
        'interpreter','latex','fontsize',18)
h=legend(h, 'correlated sources','random sources',...
    'N_{Helgason}/N_{stack}','(N_{Helgason}/N_{stack})^2');
set(h,'fontsize',8)

end
savename=strcat(savedir,dt.name,'_PS_ratio');
print(savename,'-dpng');%close


figure
h(1) = semilogx(l, ratiodat(1).corr,'r-','linewidth',3);hold on
h(2) = semilogx(l, ratiodat(1).rand,'b-','linewidth',3);hold on

for im=10:13

m_min = stackmapdat(ifield).m_min_arr(im);
m_max = stackmapdat(ifield).m_max_arr(im);
N_helg = IGLcounts_helgason(inst,(m_min+m_max)/2) * 4;
N_stack = countg_arr(im);
h1 = hline((N_helg/N_stack),'k--');
h2 = hline((N_helg/N_stack).^2,'k:');
h1.LineWidth = 2;
h2.LineWidth = 2;
end
xlim([1e2,2e5]);
title('all','fontsize',15);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$C_\ell(\rm all)/C_\ell(stacked)$',...
        'interpreter','latex','fontsize',18)
h=legend(h, 'correlated sources','random sources');
set(h,'fontsize',12)

savename=strcat(savedir,dt.name,'_PS_ratio_all');
print(savename,'-dpng');%close
%}
%%
flight=40030;
inst=1;
ifield=8;
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

srcmapdir = strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
savedir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/ihl_pl/TM',...
    num2str(inst),'/');

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
countg_arr = stackmapdat(ifield).count_stackg_arr;

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/excessdat',loaddir),'excessdat');

loaddir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/corr_sim/TM',...
    num2str(inst),'/');
load(sprintf('%smaps/ratiodat',loaddir),'ratiodat');

pltsavedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

pixscale=7;
[~,~,~,~,binl]=get_angular_spec(randn(1024),randn(1024),pixscale);

masktot = ones(1024);
srcmaptot = zeros(1024);
ihlmaptot = zeros(1024);
Nr_arr = [];
for im=10:13
    figure
    m_min = excessdat.m_min_arr(im);
    m_max = excessdat.m_max_arr(im);
        
    srcmap = fits_read(strcat(srcmapdir,dt.name,'_srcmapg',...
         num2str(m_min),'_',num2str(m_max),'ps.fits'));
    ihlmap = fits_read(strcat(savedir,dt.name,'_ihlmap_pl',...
    num2str(m_min),'_',num2str(m_max),'.fits'));
    ihlmap(find(ihlmap==Inf)) = 0;
        
    mask = make_mask_ps(flight,inst,ifield,1,m_min,m_max);
    mkk =get_mkk_sim(mask,pixscale,binl,10,numel(binl),1,ones(1024),0,NaN);

    srcmaptot = srcmaptot + srcmap;
    ihlmaptot = ihlmaptot + ihlmap;
    masktot = masktot.*mask;
    
    mmap = (srcmap+ihlmap).*mask;
    [Clihl0,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
    
    mmap = (srcmap).*mask;
    [Clpsf0,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
    
    
    N_helg = IGLcounts_helgason(inst,(m_min+m_max)/2) * 4;
    N_stack = countg_arr(im);
    Nr = N_helg./N_stack;
    Nr_arr = [Nr_arr Nr];
    
    Clihl = Clihl0.*ratiodat(im).corr;
    Clpsf = Clpsf0.*ratiodat(im).corr;
    
    h(1)=loglog(l,l.*(l+1).*Clihl0./2./pi,'b-','linewidth',2);hold on
    h(2)=loglog(l,l.*(l+1).*Clihl./2./pi,'ro-','linewidth',2);
    h(3)=loglog(l,l.*(l+1).*Clpsf./2./pi,'mo-','linewidth',2);
    h(4)=loglog(l,l.*(l+1).*Clihl0.*Nr./2./pi,'k-','linewidth',2);
    loglog(l,l.*(l+1).*Clihl0.*Nr.^2./2./pi,'k-','linewidth',2);

    xlim([1e2,2e5]);
    ylim([1e-5,1e0]);
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',20);
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
            'interpreter','latex','fontsize',18)
    leg=legend(h,{'PSF+IHL before correction','PSF+IHL sim correction',...
        'PSF sim correction','PSF+IHL correction bounds(N & N^2)'});
    set(leg,'fontsize',12);
    legend boxoff
    
    savename=strcat(pltsavedir,dt.name,'_PS',num2str(m_min),'_',num2str(m_max));
    print(savename,'-dpng');%close
end


mkk =get_mkk_sim(masktot,pixscale,binl,10,numel(binl),1,ones(1024),0,NaN);
mmap = (srcmaptot+ihlmaptot).*masktot;
[Clihl0,~,~,l]=get_Cl(mmap,masktot,mkk,pixscale,ones(1024));

mmap = (srcmap).*mask;
[Clpsf0,~,~,l]=get_Cl(mmap,masktot,mkk,pixscale,ones(1024));
Nrmin = min(Nr_arr);
Nrmax = max(Nr_arr);
Clihl = Clihl0.*ratiodat(1).corr;
Clpsf = Clpsf0.*ratiodat(1).corr;

h(1)=loglog(l,l.*(l+1).*Clihl0./2./pi,'b-','linewidth',2);hold on
h(2)=loglog(l,l.*(l+1).*Clihl./2./pi,'ro-','linewidth',2);
h(3)=loglog(l,l.*(l+1).*Clpsf./2./pi,'mo-','linewidth',2);
h(4)=loglog(l,l.*(l+1).*Clihl0.*Nrmin./2./pi,'k-','linewidth',2);
loglog(l,l.*(l+1).*Clihl0.*Nrmax.^2./2./pi,'k-','linewidth',2);
xlim([1e2,2e5]);
ylim([1e-4,1e1]);
title('total (16<m<20)','fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
leg=legend(h,{'PSF+IHL before correction','PSF+IHL sim correction',...
    'PSF sim correction','PSF+IHL correction bounds(N & N^2)'});
set(leg,'fontsize',12);
legend boxoff

savename=strcat(pltsavedir,dt.name,'_PStot');
print(savename,'-dpng');%close
