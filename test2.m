flight=40030;
inst=1;
masklim = false;
mypaths = get_paths(flight);
dx = 1200;
nbins = 25;
%% get the HSC clustering
[psf_arr,~,~] =  PSF_stacked_profile(flight,inst,8);
figure
for im=3:4
    excb_all = zeros([12,numel(psf_arr)]);
    excb_err_all = zeros([12,numel(psf_arr)]);
    for hsc_idx=0:11
        [name,~] = HSC_fields_info(hsc_idx);
        loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
        if masklim
            load(sprintf('%s/hsc/stackdathsc_%s_masklim',...
                loaddir,name),'stackdathsc');  
            load(sprintf('%s/hsc/stackdathsc_%s_masklim_bk',...
                loaddir,name),'stackdathscbk');  
        else
            load(sprintf('%s/hsc/stackdathsc_%s',...
                loaddir,name),'stackdathsc');
            load(sprintf('%s/hsc/stackdathsc_%s_bk',...
                loaddir,name),'stackdathscbk');
        end
        
        stackdat = stackdathsc(im).stackdat;
        m_min = stackdat.m_min;
        m_max = stackdat.m_max;
        r_arr = stackdat.r_arr;
        bkallcb_arr = ones([numel(stackdathscbk),numel(r_arr)]).*nan;
        for iter=1:numel(stackdathscbk)
            stackdatbk = stackdathscbk(iter).stackdat(im);
            sp  = find(stackdatbk.hitmapg_arr~=0);
            bkallcb_arr(iter,sp) = stackdatbk.profcbg_arr(sp);
        end
        bkavgcb_arr = nanmean(bkallcb_arr);
        sp = find(bkavgcb_arr==bkavgcb_arr);
        bkavgcb_arr = spline(r_arr(sp),bkavgcb_arr(sp),r_arr);
        bkavgcb_err = nanstd(bkallcb_arr);
        sp = find(bkavgcb_err==bkavgcb_err);
        bkavgcb_err = spline(r_arr(sp),bkavgcb_err(sp),r_arr);

        profcbg = stackdat.all.profcbg - bkavgcb_arr;
        profcbg_err = sqrt(stackdat.errjack.profcbg.^2 + bkavgcb_err.^2);
        excb_all(hsc_idx+1,:) = profcbg - psf_arr.*profcbg(1);
        excb_err_all(hsc_idx+1,:) = profcbg_err;        
    end
    excb_all = sum(excb_all./(excb_err_all.^2))./sum(1./(excb_err_all.^2));
    excb_err_all = sqrt(1./sum(1./(excb_err_all.^2)));
    h=loglog(r_arr,excb_all,'linewidth',2,...
        'DisplayName',sprintf('data: %d < m < %d',m_min,m_max));hold on
    c = polyfit(log10(r_arr(r_arr > 7)),log10(excb_all(r_arr > 7)),1);
    clusfit = 10.^polyval(c,log10(r_arr(r_arr > 7)));
    loglog(r_arr(r_arr > 7), clusfit,'--','color',get(h,'Color'),'linewidth',2,...
        'DisplayName',sprintf('fit: I = %.3fr^{%.2f}',10^c(2),c(1)));
end
xlim([4e-1,1.1e3])
ylim([3e-3,1.1])
title('HSC Clustering Excess');
ylabel('Excess I [nW/m^2/sr]', 'fontsize',15)
xlabel('r [arcsec]', 'fontsize',15)
h=legend('show','Location','northeast');
set(h,'fontsize',15)
legend boxoff

% binning clus to 15 rad bins
clushscfull = zeros(size(r_arr));
clushscfull(r_arr > 7) = clusfit;
dt=get_dark_times(flight,inst,8);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if masklim
    load(sprintf('%s/stackdat_%s_masklim',...
            loaddir,dt.name),'stackdatall');        
else
    load(sprintf('%s/stackdat_%s',...
            loaddir,dt.name),'stackdatall');
end
w = stackdatall(4).stackdat.radweight;
clushsc = zeros([1,15]);
clushsc(2:end-1) = clushscfull(7:19);
clushsc(1) = sum(clushscfull(1:6).*w(1:6))./sum(w(1:6));
clushsc(end) = sum(clushscfull(20:25).*w(20:25))./sum(w(20:25));
%% plot the example function param dependence
ifield = 8;
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if masklim
    load(sprintf('%s/stackdat_%s_masklim',...
            loaddir,dt.name),'stackdatall');        
else
    load(sprintf('%s/stackdat_%s',...
            loaddir,dt.name),'stackdatall');
end
r_arr = stackdatall(1).stackdat.r_arr;
psfwin_arr0 = stackdatall(1).stackdat.psfcb;

radmap = make_radius_map(zeros(2*dx+1),dx+1,dx+1).*0.7;
psfwin_map = spline(r_arr(1:18),psfwin_arr0(1:18),radmap);
psfwin_map(radmap > r_arr(18)) = 0;
profile = radial_prof(psfwin_map,ones(2*dx+1),dx+1,dx+1,1,nbins);
psfwin_arr = profile.prof./profile.prof(1);

im = 2;
rrmap = make_radius_map(zeros(201),101,101).*0.7;
[hsc_map,hsc_map_sub,R200,param_arr] = HSC_Wang19_prof(rrmap,im,true);
hsc_map1 = squeeze(hsc_map_sub(1,:,:));
hsc_map2 = squeeze(hsc_map_sub(2,:,:));

psfwinhsc_map1 = conv2(psfwin_map,hsc_map1,'same');
profile = radial_prof(psfwinhsc_map1,ones(2*dx+1),dx+1,dx+1,1,nbins);
psfwinhsc1_arr = profile.prof;

psfwinhsc_map2 = conv2(psfwin_map,hsc_map2,'same');
profile = radial_prof(psfwinhsc_map2,ones(2*dx+1),dx+1,dx+1,1,nbins);
psfwinhsc2_arr = profile.prof;

%%% get data and clus %%%
clushsc_norm = clushsc./stackdatall(im).stackdat.all.profcbg(1);
rsub_arr = stackdatall(im).stackdat.rsub_arr;
w = stackdatall(im).stackdat.radweight;
datfull = stackdatall(im).stackdat.excess.diffcb;
dat = zeros([1,15]);
dat(2:end-1) = datfull(7:19);
dat(1) = sum(datfull(1:6).*w(1:6))./sum(w(1:6));
dat(end) = sum(datfull(20:25).*w(20:25))./sum(w(20:25));
dat = dat./stackdatall(im).stackdat.all.profcbg(1);
dat_err = sqrt(diag(stackdatall(im).stackdat.cov.cov_matsub)')...
    ./stackdatall(im).stackdat.all.profcbg(1);
%%%%% varying Ie %%%%%
figure
setwinsize(gcf,1500,300)
subplot(1,4,2)
loglog(rsub_arr,dat,'k.','markersize',10);hold on
loglog(rsub_arr,-dat,'ko','markersize',5);
loglog(rsub_arr,clushsc_norm,'m--','linewidth',2);
h=legend({'CIBER excess (+)', 'CIBER excess (-)', 'HSC clustering'},...
    'Location','southwest');
set(h,'fontsize',10)
legend boxoff
errorbar(rsub_arr, dat, dat_err,'k.','markersize',10);
errorbar(rsub_arr, -dat, dat_err,'ko','markersize',5);
for Ie=[-8,-9,-10]
    hsc_map2i = sersic(rrmap./R200,Ie,param_arr(2,2),param_arr(2,3)+param_arr(1,3));
    psfwinhsc2_mapi = conv2(psfwin_map, hsc_map2i,'same');
    profile = radial_prof(psfwinhsc2_mapi,ones(2*dx+1),dx+1,dx+1,1,nbins);
    psfwinhsc2_arri = profile.prof;
    norm = psfwinhsc1_arr(1) + psfwinhsc2_arri(1);
    e = (psfwinhsc1_arr+psfwinhsc2_arri)./norm - psfwin_arr;
    subplot(1,4,1)
    h = loglog(r_arr,(psfwinhsc1_arr+psfwinhsc2_arri)./norm,'linewidth',2,...
        'DisplayName',sprintf('Ie = %d',Ie));hold on
    subplot(1,4,2)
    loglog(r_arr,e,'linewidth',2);hold on
    subplot(1,4,3)
    loglog(r_arr,psfwinhsc2_arri./norm,'color',get(h,'Color'),'linewidth',2);hold on
    subplot(1,4,4)
    loglog(r_arr,psfwinhsc1_arr./norm,'color',get(h,'Color'),'linewidth',2);hold on
end

norm = psfwinhsc1_arr(1) + psfwinhsc2_arr(1);
subplot(1,4,1)
loglog(r_arr,(psfwinhsc1_arr+psfwinhsc2_arr)./norm,'k','linewidth',1,...
    'DisplayName',sprintf('Ie = %.2f',param_arr(2,1)));
loglog(r_arr,psfwin_arr,'k--','linewidth',1,'DisplayName','PSF');hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
xlabel('r [arcsec]', 'fontsize',15)
ylabel('I / I(r=0)', 'fontsize',15)
title('Profile');
h=legend('show','Location','southwest');
set(h,'fontsize',10)
legend boxoff

subplot(1,4,2)
e = (psfwinhsc1_arr+psfwinhsc2_arr)./norm - psfwin_arr;
loglog(r_arr,e,'k','linewidth',1);hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
title('Excess');
xlabel('r [arcsec]', 'fontsize',15)

subplot(1,4,3)
loglog(r_arr,psfwinhsc2_arr./norm,'k','linewidth',1);hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
title('Outer Sersic Component');
xlabel('r [arcsec]', 'fontsize',15)

subplot(1,4,4)
loglog(r_arr,psfwinhsc1_arr./norm,'k','linewidth',1);hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
title('Inner Sersic Component');
xlabel('r [arcsec]', 'fontsize',15)

%%%%% varying n %%%%%
figure
setwinsize(gcf,1500,300)
subplot(1,4,2)
loglog(rsub_arr,dat,'k.','markersize',10);hold on
loglog(rsub_arr,-dat,'ko','markersize',5);
loglog(rsub_arr,clushsc_norm,'m--','linewidth',2);
h=legend({'CIBER excess (+)', 'CIBER excess (-)', 'HSC clustering'},...
    'Location','southwest');
set(h,'fontsize',10)
legend boxoff
errorbar(rsub_arr, dat, dat_err,'k.','markersize',10);
errorbar(rsub_arr, -dat, dat_err,'ko','markersize',5);
for n=[2,2.5,3]
    hsc_map2i = sersic(rrmap./R200,param_arr(2,1),n,param_arr(2,3)+param_arr(1,3));
    psfwinhsc2_mapi = conv2(psfwin_map, hsc_map2i,'same');
    profile = radial_prof(psfwinhsc2_mapi,ones(2*dx+1),dx+1,dx+1,1,nbins);
    psfwinhsc2_arri = profile.prof;
    norm = psfwinhsc1_arr(1) + psfwinhsc2_arri(1);
    e = (psfwinhsc1_arr+psfwinhsc2_arri)./norm - psfwin_arr;
    subplot(1,4,1)
    h = loglog(r_arr,(psfwinhsc1_arr+psfwinhsc2_arri)./norm,'linewidth',2,...
        'DisplayName',sprintf('n = %.1f',n));hold on
    subplot(1,4,2)
    loglog(r_arr,e,'linewidth',2);hold on
    subplot(1,4,3)
    loglog(r_arr,psfwinhsc2_arri./norm,'color',get(h,'Color'),'linewidth',2);hold on
    subplot(1,4,4)
    loglog(r_arr,psfwinhsc1_arr./norm,'color',get(h,'Color'),'linewidth',2);hold on
end

norm = psfwinhsc1_arr(1) + psfwinhsc2_arr(1);
subplot(1,4,1)
loglog(r_arr,(psfwinhsc1_arr+psfwinhsc2_arr)./norm,'k','linewidth',1,...
    'DisplayName',sprintf('n = %.2f',param_arr(2,2)));
loglog(r_arr,psfwin_arr,'k--','linewidth',1,'DisplayName','PSF');hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
xlabel('r [arcsec]', 'fontsize',15)
ylabel('I / I(r=0)', 'fontsize',15)
title('Profile');
h=legend('show','Location','southwest');
set(h,'fontsize',10)
legend boxoff

subplot(1,4,2)
e = (psfwinhsc1_arr+psfwinhsc2_arr)./norm - psfwin_arr;
loglog(r_arr,e,'k','linewidth',1);hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
title('Excess');
xlabel('r [arcsec]', 'fontsize',15)

subplot(1,4,3)
loglog(r_arr,psfwinhsc2_arr./norm,'k','linewidth',1);hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
title('Outer Sersic Component');
xlabel('r [arcsec]', 'fontsize',15)

subplot(1,4,4)
loglog(r_arr,psfwinhsc1_arr./norm,'k','linewidth',1);hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
title('Inner Sersic Component');
xlabel('r [arcsec]', 'fontsize',15)

%%%%% varying xe %%%%%
figure
setwinsize(gcf,1500,300)
subplot(1,4,2)
loglog(rsub_arr,dat,'k.','markersize',10);hold on
loglog(rsub_arr,-dat,'ko','markersize',5);
loglog(rsub_arr,clushsc_norm,'m--','linewidth',2);
h=legend({'CIBER excess (+)', 'CIBER excess (-)', 'HSC clustering'},...
    'Location','southwest');
set(h,'fontsize',10)
legend boxoff
errorbar(rsub_arr, dat, dat_err,'k.','markersize',10);
errorbar(rsub_arr, -dat, dat_err,'ko','markersize',5);
for xe=[0.01,0.03,0.05]
    hsc_map2i = sersic(rrmap./R200,param_arr(2,1),param_arr(2,2),param_arr(1,3)+xe);
    psfwinhsc2_mapi = conv2(psfwin_map, hsc_map2i,'same');
    profile = radial_prof(psfwinhsc2_mapi,ones(2*dx+1),dx+1,dx+1,1,nbins);
    psfwinhsc2_arri = profile.prof;
    norm = psfwinhsc1_arr(1) + psfwinhsc2_arri(1);
    e = (psfwinhsc1_arr+psfwinhsc2_arri)./norm - psfwin_arr;
    subplot(1,4,1)
    h = loglog(r_arr,(psfwinhsc1_arr+psfwinhsc2_arri)./norm,'linewidth',2,...
        'DisplayName',sprintf('xe = %.2f',xe));hold on
    subplot(1,4,2)
    loglog(r_arr,e,'linewidth',2);hold on
    subplot(1,4,3)
    loglog(r_arr,psfwinhsc2_arri./norm,'color',get(h,'Color'),'linewidth',2);hold on
    subplot(1,4,4)
    loglog(r_arr,psfwinhsc1_arr./norm,'color',get(h,'Color'),'linewidth',2);hold on
end

norm = psfwinhsc1_arr(1) + psfwinhsc2_arr(1);
subplot(1,4,1)
loglog(r_arr,(psfwinhsc1_arr+psfwinhsc2_arr)./norm,'k','linewidth',1,...
    'DisplayName',sprintf('xe = %.2f',param_arr(2,3)));
loglog(r_arr,psfwin_arr,'k--','linewidth',1,'DisplayName','PSF');hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
xlabel('r [arcsec]', 'fontsize',15)
ylabel('I / I(r=0)', 'fontsize',15)
title('Profile');
h=legend('show','Location','southwest');
set(h,'fontsize',10)
legend boxoff

subplot(1,4,2)
e = (psfwinhsc1_arr+psfwinhsc2_arr)./norm - psfwin_arr;
loglog(r_arr,e,'k','linewidth',1);hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
title('Excess');
xlabel('r [arcsec]', 'fontsize',15)

subplot(1,4,3)
loglog(r_arr,psfwinhsc2_arr./norm,'k','linewidth',1);hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
title('Outer Sersic Component');
xlabel('r [arcsec]', 'fontsize',15)

subplot(1,4,4)
loglog(r_arr,psfwinhsc1_arr./norm,'k','linewidth',1);hold on
xlim([4e-1,700])
ylim([1e-5,1.1])
title('Inner Sersic Component');
xlabel('r [arcsec]', 'fontsize',15)
%% fitting
ifield = 8;
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if masklim
    load(sprintf('%s/stackdat_%s_masklim',...
            loaddir,dt.name),'stackdatall');        
else
    load(sprintf('%s/stackdat_%s',...
            loaddir,dt.name),'stackdatall');
end
r_arr = stackdatall(1).stackdat.r_arr;
psfwin_arr0 = stackdatall(1).stackdat.psfcb;

radmap = make_radius_map(zeros(2*dx+1),dx+1,dx+1).*0.7;
psfwin_map = spline(r_arr(1:18),psfwin_arr0(1:18),radmap);
psfwin_map(radmap > r_arr(18)) = 0;
profile = radial_prof(psfwin_map,ones(2*dx+1),dx+1,dx+1,1,nbins);
psfwin_arr = profile.prof./profile.prof(1);

im = 2;
rrmap = make_radius_map(zeros(201),101,101).*0.7;
[hsc_map,hsc_map_sub,R200,param_arr] = HSC_Wang19_prof(rrmap,im,true);
hsc_map1 = squeeze(hsc_map_sub(1,:,:));
hsc_map2 = squeeze(hsc_map_sub(2,:,:));

psfwinhsc_map1 = conv2(psfwin_map,hsc_map1,'same');
profile = radial_prof(psfwinhsc_map1,ones(2*dx+1),dx+1,dx+1,1,nbins);
psfwinhsc1_arr = profile.prof;

psfwinhsc_map2 = conv2(psfwin_map,hsc_map2,'same');
profile = radial_prof(psfwinhsc_map2,ones(2*dx+1),dx+1,dx+1,1,nbins);
psfwinhsc2_arr = profile.prof;

%%% get data and clus %%%
rsub_arr = stackdatall(im).stackdat.rsub_arr;
clushsc_norm = clushsc./stackdatall(im).stackdat.all.profcbg(1);
w = stackdatall(im).stackdat.radweight;
datfull = stackdatall(im).stackdat.excess.diffcb;
dat = zeros([1,15]);
dat(2:end-1) = datfull(7:19);
dat(1) = sum(datfull(1:6).*w(1:6))./sum(w(1:6));
dat(end) = sum(datfull(20:25).*w(20:25))./sum(w(20:25));
dat = dat./stackdatall(im).stackdat.all.profcbg(1);
dat_err = sqrt(diag(stackdatall(im).stackdat.cov.cov_matsub)')...
    ./stackdatall(im).stackdat.all.profcbg(1);
cov_inv = stackdatall(im).stackdat.cov.cov_invsub;
cov_inv = cov_inv.*(stackdatall(im).stackdat.all.profcbg(1)^2);
figure
for xe=[0.01,0.02,0.03,0.05]
    hsc_map2i = sersic(rrmap./R200,param_arr(2,1),param_arr(2,2),param_arr(1,3)+xe);
    psfwinhsc2_mapi = conv2(psfwin_map, hsc_map2i,'same');
    profile = radial_prof(psfwinhsc2_mapi,ones(2*dx+1),dx+1,dx+1,1,nbins);
    psfwinhsc2_arri = profile.prof;
    norm = psfwinhsc1_arr(1) + psfwinhsc2_arri(1);
    efull = (psfwinhsc1_arr+psfwinhsc2_arri)./norm - psfwin_arr;
    e = zeros([1,15]);
    e(2:end-1) = efull(7:19);
    e(1) = sum(efull(1:6).*w(1:6))./sum(w(1:6));
    e(end) = sum(efull(20:25).*w(20:25))./sum(w(20:25));
    for clus_scale = [1,2,3]
        m = clus_scale.*clushsc_norm + e;
        diff = dat - m;
        chi2 = diff*cov_inv*diff';
        chi2sub = diff(1:14)*cov_inv(1:14,1:14)*diff(1:14)';
        fprintf('xe = %.3f, clus_scale = %d,chi2 = %.2e,chi2sub = %.2e\n',...
            xe,clus_scale,chi2,chi2sub); 
        if chi2sub > 200
            continue
        end
        loglog(rsub_arr,m,'DisplayName',...
        sprintf('(xe,clus scale) = (%.2f, %.1f) \n chi^2 = %.1f (%.1f)',...
        xe,clus_scale,chi2,chi2sub));hold on
    end
end
h=legend('show','Location','southwest');
set(h,'fontsize',10)
legend boxoff
loglog(rsub_arr,dat,'k.','markersize',10);hold on
loglog(rsub_arr,-dat,'ko','markersize',5);
errorbar(rsub_arr, dat, dat_err,'k.','markersize',10);
errorbar(rsub_arr, -dat, dat_err,'ko','markersize',5);
xlabel('r [arcsec]', 'fontsize',15)
ylabel('I / I(r=0)', 'fontsize',15)
xlim([4e-1,700])
ylim([1e-5,1.1])