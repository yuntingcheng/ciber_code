function plot_ihl_stack_hist(flight,inst,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addOptional('rmin',nan,@isnumeric);
  
  p.parse(flight,inst,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  rmin     = p.Results.rmin;
  
  clear p varargin;

%%
mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

%%
%%% plot the stacking profile %%%
for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if rmin==2
    load(sprintf('%s/%s_ihlprofdat_hist_rmin2',loaddir,dt.name),'ihlprofdat');
elseif isnan(rmin)
    load(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');
end

figure
setwinsize(gcf,1200,750)
for im=1:4
    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;
    counts = ihlprofdat.data(im).counts;
    countg = ihlprofdat.data(im).countg;
    r_arr = ihlprofdat.r_arr;

    profscb = ihlprofdat.data(im).profscb;
    profscb_err = ihlprofdat.data(im).profscb_err;
    profsps = ihlprofdat.data(im).profsps;
    profsps_err = ihlprofdat.data(im).profsps_err;
    profgcb = ihlprofdat.data(im).profgcb;
    profgcb_err = ihlprofdat.data(im).profgcb_err;
    profgps = ihlprofdat.data(im).profgps;
    profgps_err = ihlprofdat.data(im).profgps_err;

    profscbbk = ihlprofdat.bk(im).profscb;
    profscb_errbk = ihlprofdat.bk(im).profscb_err;
    profspsbk = ihlprofdat.bk(im).profsps;
    profsps_errbk = ihlprofdat.bk(im).profsps_err;

    subplot(3,4,im)
    semilogx(r_arr.*0.98,profscb,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgcb,'b','markersize',10);
    plot(r_arr.*0.98,profsps,'m','markersize',10);
    plot(r_arr.*1.02,profgps,'c','markersize',10);
    h=legend({'CIBER stars','CIBER galaxies',...
        'PanSTARRS stars','PanSTARRS galaxies'},'Location','southwest');
    set(h,'fontsize',7)
    legend boxoff
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    errorbar(r_arr.*0.98,profsps,profsps_err,'m.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'c.','markersize',10);
    xlim([4e-1,50])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
    title({strcat(num2str(m_min),'<m<',num2str(m_max)),...
        strcat(num2str(counts),' stars, ',num2str(countg), ' galaxies')},...
        'fontsize',15)
            
    subplot(3,4,im + 4)
    semilogx(r_arr.*0.98,profscb,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgcb,'b','markersize',10);
    plot(r_arr,profscbbk,'k','markersize',10);
    h=legend({'CIBER stars','CIBER galaxies','background'},...
        'Location','northeast');
    set(h,'fontsize',7)
    legend boxoff
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    errorbar(r_arr,profscbbk,profscb_errbk,'k.','markersize',10);
    xlim([4e-1,1e3])
    ylim([-1,10])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)

    subplot(3,4,im + 8)
    semilogx(r_arr.*0.98,profsps,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgps,'b','markersize',10);
    plot(r_arr,profspsbk,'k','markersize',10);
    h=legend({'PanSTARRS stars','PanSTARRS galaxies'},...
        'Location','southeast');
    set(h,'fontsize',7)
    legend boxoff
    errorbar(r_arr.*0.98,profsps,profsps_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'b.','markersize',10);
    errorbar(r_arr,profspsbk,profsps_errbk,'k.','markersize',10);
    xlim([4e-1,1e3])
    ylim([-0.05,0.1])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end
suptitle(dt.name);

if rmin==2
    savename=strcat(pltsavedir,dt.name,'_stackprof_hist_rmin2');
elseif isnan(rmin)
    savename=strcat(pltsavedir,dt.name,'_stackprof_hist');
end

print(savename,'-dpng');%close
end

%%
%%% plot the excess profile in each field %%%

for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if rmin==2
    load(sprintf('%s/%s_ihlprofdat_hist_rmin2',loaddir,dt.name),'ihlprofdat');
elseif isnan(rmin)
    load(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');
end

figure
setwinsize(gcf,1400,360)
for im=1:4
    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;
    counts = ihlprofdat.data(im).counts;
    countg = ihlprofdat.data(im).countg;
    r_arr = ihlprofdat.r_arr;

    ecb = ihlprofdat.excess(im).diffcb;
    ecb_err = ihlprofdat.excess(im).diffcb_err1;
    eps = ihlprofdat.excess(im).diffps;
    eps_err = ihlprofdat.excess(im).diffps_err1;
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;
    
    subplot(1,4,im)
    loglog(r_arr.*0.98,ecb,'r.','markersize',10);hold on
    loglog(r_arr.*1.02,eps,'b.','markersize',10);hold on
    loglog(r_arr,ihl,'k.','markersize',10);hold on
    
    sp = find(r_arr > 100);
    v = ihlprofdat.excess(im).diff100;
    e = ihlprofdat.excess(im).diff_err100;
    if v>=0 & v-e >=0
        fill([r_arr(sp),flip(r_arr(sp))],...
            [(v+e) * ones(size(sp)),(v-e) * ones(size(sp))],...
                [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
        plot(r_arr(sp),v*ones(size(sp)),'color',[0.2,0.4,0.2],'linewidth',2);
    elseif v>=0 & v-e<0
        fill([r_arr(sp),flip(r_arr(sp))],...
            [(v+e) * ones(size(sp)),1e-8 * ones(size(sp))],...
                [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
        plot(r_arr(sp),v*ones(size(sp)),'color',[0.2,0.4,0.2],'linewidth',2);
    elseif v<0 & v+e>0
        fill([r_arr(sp),flip(r_arr(sp))],...
            [-(v-e) * ones(size(sp)),1e-8 * ones(size(sp))],...
                [0.6,0.4,0.2],'facealpha',0.5,'EdgeColor','none');
        plot(r_arr(sp),-v*ones(size(sp)),'--','color',[0.6,0.4,0.2],'linewidth',2);
    elseif v<0 & v+e<=0
        fill([r_arr(sp),flip(r_arr(sp))],...
            [-(v-e) * ones(size(sp)),-(v+e) * ones(size(sp))],...
                [0.6,0.4,0.2],'facealpha',0.5,'EdgeColor','none');
        plot(r_arr(sp),-v*ones(size(sp)),'--','color',[0.6,0.4,0.2],'linewidth',2);
    end

    h=legend({'CIBER Excess','PanSTARRS Excess', ...
        'Excess','Excess > 100 arcsec'},'Location','northeast');
    set(h,'fontsize',10)
    legend boxoff

    eps_err_low = eps_err;
    for ii=1:numel(eps)
        if eps(ii)>0 & eps_err(ii) > eps(ii)
            eps_err_low(ii) = eps(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*1.02,eps,eps_err_low,eps_err,'b.','markersize',10);
    loglog(r_arr.*1.02,-eps,'bo','markersize',5);hold on
    eps_err_low = eps_err;
    for ii=1:numel(eps)
        if eps(ii)<0 & eps_err(ii) > abs(eps(ii))
            eps_err_low(ii) = -eps(ii) - 1e-8;
        end
    end    
    errorbar(r_arr.*1.02,-eps,eps_err_low,eps_err,'bo','markersize',5);
    
    ecb_err_low = ecb_err;
    for ii=1:numel(eps)
        if ecb(ii)>0 & ecb_err(ii) > ecb(ii)
            ecb_err_low(ii) = ecb(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*0.98,ecb,ecb_err_low,ecb_err,'r.','markersize',10);
    loglog(r_arr.*0.98,-ecb,'ro','markersize',5);
    ecb_err_low = ecb_err;
    for ii=1:numel(ecb)
        if ecb(ii)<0 & ecb_err(ii) > abs(ecb(ii))
            ecb_err_low(ii) = -ecb(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*0.98,-ecb,ecb_err_low,ecb_err,'ro','markersize',5);
    
    ihl_err_low = ihl_err;
    for ii=1:numel(ihl)
        if ihl(ii)>0 & ihl_err(ii) > ihl(ii)
            ihl_err_low(ii) = ihl(ii) - 1e-8;
        end
    end
    errorbar(r_arr,ihl,ihl_err_low,ihl_err,'k.','markersize',10);
    loglog(r_arr,-ihl,'ko','markersize',5);
    
    ihl_err_low = ihl_err;
    for ii=1:numel(ecb)
        if ihl(ii)<0 & ihl_err(ii) > abs(ihl(ii))
            ihl_err_low(ii) = -ihl(ii) - 1e-8;
        end
    end    
    errorbar(r_arr,-ihl,ihl_err_low,ihl_err,'ko','markersize',5);
    xlim([4e-1,1e3])
    ylim([1e-4,1e4])
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
      ' (',num2str(counts),' stars, ',num2str(countg), ' galaxies)'),'fontsize',10)
    xlabel('arcsec', 'fontsize',10)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end
suptitle(dt.name);
if rmin==2
    savename=strcat(pltsavedir,dt.name,'_excessprof_hist_rmin2');
elseif isnan(rmin)
    savename=strcat(pltsavedir,dt.name,'_excessprof_hist');
end
print(savename,'-dpng');%close

end
%%
%%% plot the excess profile >100 arcsec in each field %%%

figure
setwinsize(gcf,1400,300)
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
for im=1:4
    subplot(1,4,im)
    name = {};
    for ifield = 4:8
        dt=get_dark_times(flight,inst,ifield);
        name{end+1} = dt.name;
        if rmin==2
            load(sprintf('%s/%s_ihlprofdat_hist_rmin2',loaddir,dt.name),...
                'ihlprofdat');
        elseif isnan(rmin)
            load(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');
        end

        m_min = ihlprofdat.data(im).m_min;
        m_max = ihlprofdat.data(im).m_max;
        r_arr = ihlprofdat.r_arr;

        sp = find(r_arr > 100);
        v = ihlprofdat.excess(im).diff100;
        e = ihlprofdat.excess(im).diff_err100;
        errorbar(ifield-3,v,e,'k.','markersize',10);hold on
    end
    plot([0,6],[0,0],'b--')
    ylabel('mean excess > 100 arcsec [nW/m^2/sr]', 'fontsize',10);
    xlim([0.5,5.5])
    set(gca,'xTick',1:5);
    set(gca, 'xTickLabels', name); 
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);

end

if rmin==2
    savename=strcat(pltsavedir,'excess_fields_rmin2');
elseif isnan(rmin)
    savename=strcat(pltsavedir,'excess_fields');
end
print(savename,'-dpng');%close
%%
%%% get excess profile in each field%%%

for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if rmin==2
    load(sprintf('%s/%s_ihlprofdat_hist_rmin2',loaddir,dt.name),'ihlprofdat');
elseif isnan(rmin)
    load(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');
end
for im=1:4
    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;
    counts = ihlprofdat.data(im).counts;
    countg = ihlprofdat.data(im).countg;
    r_arr = ihlprofdat.r_arr;

    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;

    excess(ifield).count(im).counts = counts;
    excess(ifield).count(im).countg = countg;
    excess(ifield).ihl(im).ihl = ihl;
    excess(ifield).ihl(im).err = ihl_err;
    excess(ifield).ihl(im).ihl100 = ihlprofdat.excess(im).diff100;
    excess(ifield).ihl(im).err100 = ihlprofdat.excess(im).diff_err100;

end
end

%%

figure
setwinsize(gcf,1500,600)

for im=1:4
ihldat = zeros([5,numel(ihl)]);
errdat = zeros([5,numel(ihl)]);

for ifield = 4:8
    dt=get_dark_times(flight,inst,ifield);

    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;

    ihl = excess(ifield).ihl(im).ihl;
    ihl_err = excess(ifield).ihl(im).err;
    ihldat(ifield-3,:) = ihl;
    errdat(ifield-3,:) = ihl_err;
    
    off = 1 + (ifield - 6).*0.02;
    subplot(2,4,im)
    loglog(r_arr.*off,ihl,'.','color',...
        get_color(ifield-3),'markersize',10,'DisplayName',dt.name);hold on
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);

    
end
h=legend('show','Location','southwest');
set(h,'fontsize',10)
legend boxoff
ihlmid = median(ihldat);
errmid = median(errdat);

for ifield = 4:8
    ihl = excess(ifield).ihl(im).ihl;
    ihl_err = excess(ifield).ihl(im).err; 
    off = 1 + (ifield - 6).*0.02;
    subplot(2,4,im)
    loglog(r_arr.*off,ihl,'.','color',...
        get_color(ifield-3),'markersize',10);hold on
    errorbar(r_arr.*off,ihl,ihl_err,'.','color',...
        get_color(ifield-3),'markersize',10);
    loglog(r_arr.*off,-ihl,'o','color',...
        get_color(ifield-3),'markersize',5);hold on
    errorbar(r_arr.*off,-ihl,ihl_err,'o','color',...
        get_color(ifield-3),'markersize',5);
    xlim([4e-1,1e3])
    ylim([1e-4,1e3])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)

end

for ifield = 4:8

    ihl = excess(ifield).ihl(im).ihl;
    ihl_err = excess(ifield).ihl(im).err;

    off = 1 + (ifield - 6).*0.02;
    subplot(2,4,im + 4)
    semilogx(r_arr.*off,(ihl - ihlmid)./errmid,'.','color',...
        get_color(ifield-3),'markersize',10);hold on
    errorbar(r_arr.*off,(ihl - ihlmid)./errmid,ihl_err./errmid,'.','color',...
        get_color(ifield-3),'markersize',10);
    xlim([4e-1,1e3])
    xlabel('arcsec', 'fontsize',15)
    ylabel('scaled excess', 'fontsize',15)

end
end

if rmin==2
    savename=strcat(pltsavedir,'excess_all_hist_rmin2');
elseif isnan(rmin)
    savename=strcat(pltsavedir,'excess_all_hist');
end
print(savename,'-dpng');%close
%%
figure
setwinsize(gcf,1400,300)

for im=1:4
ihldat = zeros([5,numel(ihl)]);
errdat = zeros([5,numel(ihl)]);
ihldat100 = zeros([5,1]);
errdat100 = zeros([5,1]);
m_min = ihlprofdat.data(im).m_min;
m_max = ihlprofdat.data(im).m_max;
for ifield = 4:8

    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;

    ihl = excess(ifield).ihl(im).ihl;
    ihl_err = excess(ifield).ihl(im).err;
    ihldat(ifield-3,:) = ihl;
    errdat(ifield-3,:) = ihl_err;
    ihldat100(ifield-3,:) = excess(ifield).ihl(im).ihl100;
    errdat100(ifield-3,:) = excess(ifield).ihl(im).err100;
end
ihltot = sum(ihldat./errdat.^2)./sum(1./errdat.^2);
ihltot_err = sqrt(1./(sum(1./errdat.^2)));
excess(1).ihl.ihl = ihltot;
excess(1).ihl.err = ihltot_err;

sp = find(r_arr > 100);
v = sum(ihldat100./errdat100.^2)./sum(1./errdat100.^2);
e = sqrt(1./sum(1./errdat100.^2));

subplot(1,4,im)
loglog(r_arr,ihltot,'k.','markersize',10);hold on
if v>=0 & v-e >=0
    fill([r_arr(sp),flip(r_arr(sp))],...
        [(v+e) * ones(size(sp)),(v-e) * ones(size(sp))],...
            [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
    plot(r_arr(sp),v*ones(size(sp)),'color',[0.2,0.4,0.2],'linewidth',2);
elseif v>=0 & v-e<0
    fill([r_arr(sp),flip(r_arr(sp))],...
        [(v+e) * ones(size(sp)),1e-8 * ones(size(sp))],...
            [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
    plot(r_arr(sp),v*ones(size(sp)),'color',[0.2,0.4,0.2],'linewidth',2);
elseif v<0 & v+e>0
    fill([r_arr(sp),flip(r_arr(sp))],...
        [-(v-e) * ones(size(sp)),1e-8 * ones(size(sp))],...
            [0.6,0.4,0.2],'facealpha',0.5,'EdgeColor','none');
    plot(r_arr(sp),-v*ones(size(sp)),'--','color',[0.6,0.4,0.2],'linewidth',2);
elseif v<0 & v+e<=0
    fill([r_arr(sp),flip(r_arr(sp))],...
        [-(v-e) * ones(size(sp)),-(v+e) * ones(size(sp))],...
            [0.6,0.4,0.2],'facealpha',0.5,'EdgeColor','none');
    plot(r_arr(sp),-v*ones(size(sp)),'--','color',[0.6,0.4,0.2],'linewidth',2);
end
legend({'averaged exess','excess > 100 arcsec'})
legend boxoff

ihltot_err_low = ihltot_err;
for ii=1:numel(eps)
    if ihltot(ii)>=0 & ihltot_err(ii) > ihltot(ii)
        ihltot_err_low(ii) = ihltot(ii) - 1e-8;
    end
end
errorbar(r_arr,ihltot,ihltot_err_low,ihltot_err,'k.','markersize',10);
loglog(r_arr,-ihltot,'ko','markersize',5);hold on
ihltot_err_low = ihltot_err;
for ii=1:numel(eps)
    if ihltot(ii)<0 & ihltot_err(ii) > abs(ihltot(ii))
        ihltot_err_low(ii) = -ihltot(ii) - 1e-8;
    end
end
errorbar(r_arr,-ihltot,ihltot_err_low,ihltot_err,'ko','markersize',5);

ylim([1e-4,1e3])
xlim([4e-1,1e3])
title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);
xlabel('arcsec', 'fontsize',15);
ylabel('I [nW/m^2/sr]', 'fontsize',15);

end
if rmin==2
    savename=strcat(pltsavedir,'excess_avg_hist_rmin2');
elseif isnan(rmin)
    savename=strcat(pltsavedir,'excess_avg_hist');
end
print(savename,'-dpng');%close
%% chi & chi2
%{
figure
setwinsize(gcf,1500,600)

for im=1:3
m_min = ihlprofdat.data(im).m_min;
m_max = ihlprofdat.data(im).m_max;

rmask = 0;
for ifield = 4:8
    rmask = rmask + get_mask_radius(inst,ifield,(m_min + m_max)/2);
end
rmask = rmask ./ 5;
sp1 = find(r_arr > rmask & r_arr <= 100);
sp2 = find(r_arr > 100);

chi1dat_arr = zeros([3,5]);
chi2dat_arr = zeros([3,5]);
name = {};
mu=excess(1).ihl.ihl;
for ifield = 4:8
    dt=get_dark_times(flight,inst,ifield);
    name{end+1} = dt.name;
    diff = excess(ifield).ihl.ihl;
    diff_err = excess(ifield).ihl.err;
    chi1dat_arr(1,ifield-3) = sum((diff - mu) ./ diff_err);
    chi2dat_arr(1,ifield-3) = sum((diff - mu).^2 ./ diff_err.^2);
    chi1dat_arr(2,ifield-3) = sum((diff(sp1) - mu(sp1)) ./ diff_err(sp1));
    chi2dat_arr(2,ifield-3) = sum((diff(sp1) - mu(sp1)).^2 ./ diff_err(sp1).^2);
    chi1dat_arr(3,ifield-3) = sum((diff(sp2) - mu(sp2)) ./ diff_err(sp2));
    chi2dat_arr(3,ifield-3) = sum((diff(sp2) - mu(sp2)).^2 ./ diff_err(sp2).^2);
end
legend_name = {};
legend_name{end+1} = sprintf('all (%d bins)',numel(r_arr));
legend_name{end+1} = sprintf('r(mask) < r < 100\''\'' (%d bins)',numel(sp1));
legend_name{end+1} = sprintf('r > 100\''\'' (%d bins)',numel(sp2));

subplot(2,3,im)
plot(4:8,chi1dat_arr(1,:)./numel(r_arr),'k.','markersize',20);hold on
plot(4:8,chi1dat_arr(2,:)./numel(sp1),'r.','markersize',20);
plot(4:8,chi1dat_arr(3,:)./numel(sp2),'b.','markersize',20);
legend(legend_name,'Location','southwest');
set(gca,'xTick',4:8); 
set(gca, 'xTickLabels', name); 
xlim([3.8,8.2])
title(strcat(num2str(m_min),' < m < ',num2str(m_max)),'fontsize', 15)
ylabel('\chi / N', 'fontsize', 15)
subplot(2,3,im+3)
plot(4:8,chi2dat_arr(1,:)./numel(r_arr),'k.','markersize',20);hold on
plot(4:8,chi2dat_arr(2,:)./numel(sp1),'r.','markersize',20);
plot(4:8,chi2dat_arr(3,:)./numel(sp2),'b.','markersize',20);
set(gca,'xTick',4:8);
set(gca, 'xTickLabels', name); 
xlim([3.8,8.2])
title(strcat(num2str(m_min),' < m < ',num2str(m_max)),'fontsize', 15)
ylabel('\chi^2 / N', 'fontsize', 15)
end
if rmin==2
    savename=strcat(pltsavedir,'excess_chi_hist_rmin2');
elseif isnan(rmin)
    savename=strcat(pltsavedir,'excess_chi_hist');
end
print(savename,'-dpng');%close
%}
%%
figure

for im=1:4
ihldat = zeros([5,numel(ihl)]);
errdat = zeros([5,numel(ihl)]);
off = 1 + (im - 2).*0.01;
for ifield = 4:8
    ihl = excess(ifield).ihl(im).ihl;
    ihl_err = excess(ifield).ihl(im).err;
    ihldat(ifield-3,:) = ihl;
    errdat(ifield-3,:) = ihl_err;
end
ihltot = sum(ihldat./errdat.^2)./sum(1./errdat.^2);
loglog(r_arr.*off,ihltot,'.','color',get_color(im),'markersize',10);hold on
end
legend({'16<m<17','17<m<18', '18<m<19'});
legend boxoff

for im=1:4
ihldat = zeros([5,numel(ihl)]);
errdat = zeros([5,numel(ihl)]);
off = 1 + (im - 2).*0.01;
for ifield = 4:8
    ihl = excess(ifield).ihl(im).ihl;
    ihl_err = excess(ifield).ihl(im).err;
    ihldat(ifield-3,:) = ihl;
    errdat(ifield-3,:) = ihl_err;
end
ihltot = sum(ihldat./errdat.^2)./sum(1./errdat.^2);
ihltot_err = sqrt(1./(sum(1./errdat.^2)));

errorbar(r_arr.*off,ihltot,ihltot_err,'.','color',get_color(im),'markersize',10);
loglog(r_arr.*off,-ihltot,'o','color',get_color(im),'markersize',5);hold on
errorbar(r_arr.*off,-ihltot,ihltot_err,'o','color',get_color(im),'markersize',5);
end

ylim([1e-4,1e3])
xlim([4e-1,1e3])
xlabel('arcsec', 'fontsize',15);
ylabel('I [nW/m^2/sr]', 'fontsize',15);

if rmin==2
    savename=strcat(pltsavedir,'excess_avg_all_hist_rmin2');
elseif isnan(rmin)
    savename=strcat(pltsavedir,'excess_avg_all_hist');
end
print(savename,'-dpng');%close
%%
%{
if rmin==2   
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
    load(sprintf('%s/%s_ihlprofdat_hist_rmin2',loaddir,dt.name),'ihlprofdat');
    ihlprofdat2 = ihlprofdat;
    load(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');

    figure
    setwinsize(gcf,1200,400)
    for im=1:3
        m_min = ihlprofdat.data(im).m_min;
        m_max = ihlprofdat.data(im).m_max;
        r_arr = ihlprofdat.r_arr;

        subplot(2,3,im)
        semilogx(r_arr,ihlprofdat.data(im).profgcb,'b','markersize',10);hold on
        semilogx(r_arr,ihlprofdat2.data(im).profgcb,'b--','markersize',10);
        semilogx(r_arr,ihlprofdat.data(im).profscb,'r','markersize',10);
        semilogx(r_arr,ihlprofdat2.data(im).profscb,'r--','markersize',10);
        legend({'CIBER gal, original mask','CIBER gal, aggresive mask',...
            'CIBER star, original mask','CIBER star, aggresive mask'});
        xlim([10,1e3])
        ylim([-1,4])
        xlabel('arcsec', 'fontsize',15)
        ylabel('I [nW/m^2/sr]', 'fontsize',15)
        title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15)

        subplot(2,3,im+3)
        semilogx(r_arr,ihlprofdat.data(im).profgps,'b','markersize',10);hold on
        semilogx(r_arr,ihlprofdat2.data(im).profgps,'b--','markersize',10);
        semilogx(r_arr,ihlprofdat.data(im).profsps,'r','markersize',10);
        semilogx(r_arr,ihlprofdat2.data(im).profsps,'r--','markersize',10);
        legend({'PanSTARRS gal, original mask','PanSTARRS gal, aggresive mask',...
            'PanSTARRS star, original mask','PanSTARRS star, aggresive mask'});
        xlim([10,1e3])
        ylim([-0.01,0.05])
        xlabel('arcsec', 'fontsize',15)
        ylabel('I [nW/m^2/sr]', 'fontsize',15)
    end
    suptitle(dt.name);
end
end
%}
%%
return