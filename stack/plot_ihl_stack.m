function plot_ihl_stack(flight,inst)

mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');

%%% plot the stacking profile %%%
figure
setwinsize(gcf,1200,1200)
for im=1:3
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

    subplot(3,3,im)
    semilogx(r_arr.*0.98,profscb,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgcb,'b','markersize',10);
    plot(r_arr.*0.98,profsps,'m','markersize',10);
    plot(r_arr.*1.02,profgps,'c','markersize',10);
    h=legend({'CIBER stars','CIBER galaxies','Sim stars','Sim galaxies'},...
        'Location','northeast');
    set(h,'fontsize',10)
    legend boxoff
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    errorbar(r_arr.*0.98,profsps,profsps_err,'m.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'c.','markersize',10);
    xlim([4e-1,50])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
      ' (',num2str(counts),' stars, ',num2str(countg), ' galaxies)'),'fontsize',15)
            
    subplot(3,3,im + 3)
    semilogx(r_arr.*0.98,profscb,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgcb,'b','markersize',10);
    plot(r_arr,ihlprofdat.bkprofcb,'k.','markersize',10);
    h=legend({'CIBER stars','CIBER galaxies', 'CIBER background stack'},...
        'Location','northwest');
    set(h,'fontsize',10)
    legend boxoff
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    errorbar(r_arr,ihlprofdat.bkprofcb,ihlprofdat.bkprofcb_err,'k.','markersize',10);
    xlim([4e-1,1e3])
    ylim([-1,10])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)

    subplot(3,3,im + 6)
    semilogx(r_arr.*0.98,profsps,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgps,'b','markersize',10);
    plot(r_arr,ihlprofdat.bkprofps,'k.','markersize',10);
    h=legend({'Sim stars','Sim galaxies', 'Sim background stack'},...
        'Location','northwest');
    set(h,'fontsize',10)
    legend boxoff
    errorbar(r_arr.*0.98,profsps,profsps_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'b.','markersize',10);
    errorbar(r_arr,ihlprofdat.bkprofps,ihlprofdat.bkprofps_err,'k.','markersize',10);
    xlim([4e-1,1e3])
    ylim([-0.05,0.1])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end
suptitle(dt.name);

savename=strcat(pltsavedir,dt.name,'_stackprof');
print(savename,'-dpng');%close
end
%%
for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');

figure
setwinsize(gcf,1200,360)
for im=1:3
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
    
    eps = profgps - profsps;
    eps_err = sqrt(profgps_err.^2 + profsps_err.^2);
    ecb = profgcb - profscb;
    ecb_err = sqrt(profgcb_err.^2 + profscb_err.^2);
    ihl = ecb - eps;
    ihl_err = sqrt(ecb_err.^2 + eps_err.^2);
    
    subplot(1,3,im)
    loglog(r_arr.*0.98,ecb,'r.','markersize',10);hold on
    loglog(r_arr.*1.02,eps,'b.','markersize',10);hold on
    loglog(r_arr,ihl,'k.','markersize',10);hold on
    h=legend({'CIBER gals - stars','Sim gals - stars', ...
        'CIBER excess - Sim excess'},'Location','northeast');
    set(h,'fontsize',10)
    legend boxoff

    errorbar(r_arr.*1.02,eps,eps_err,'b.','markersize',10);
    loglog(r_arr.*1.02,-eps,'bo','markersize',5);hold on
    errorbar(r_arr.*1.02,-eps,eps_err,'bo','markersize',5);
    errorbar(r_arr.*0.98,ecb,ecb_err,'r.','markersize',10);
    loglog(r_arr.*0.98,-ecb,'ro','markersize',5);hold on
    errorbar(r_arr.*0.98,-ecb,ecb_err,'ro','markersize',5);
    errorbar(r_arr,ihl,ihl_err,'k.','markersize',10);
    loglog(r_arr,-ihl,'ko','markersize',5);hold on
    errorbar(r_arr,-ihl,ihl_err,'ko','markersize',5);
    xlim([4e-1,1e3])
    ylim([1e-4,1e4])
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
      ' (',num2str(counts),' stars, ',num2str(countg), ' galaxies)'),'fontsize',15)
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end
suptitle(dt.name);
savename=strcat(pltsavedir,dt.name,'_excessprof');
print(savename,'-dpng');%close

end
%%
for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');
for im=1:3
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
    
    eps = profgps - profsps;
    eps_err = sqrt(profgps_err.^2 + profsps_err.^2);
    ecb = profgcb - profscb;
    ecb_err = sqrt(profgcb_err.^2 + profscb_err.^2);
    ihl = ecb - eps;
    ihl_err = sqrt(ecb_err.^2 + eps_err.^2);
    
    excess(ifield).count(im).counts = counts;
    excess(ifield).count(im).countg = countg;
    excess(ifield).ihl(im).ihl = ihl;
    excess(ifield).ihl(im).err = ihl_err;
end
end

figure
setwinsize(gcf,1500,600)

for im=1:3
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
    subplot(2,3,im)
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
    subplot(2,3,im)
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

    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;

    ihl = excess(ifield).ihl(im).ihl;
    ihl_err = excess(ifield).ihl(im).err;
    
    off = 1 + (ifield - 6).*0.02;
    subplot(2,3,im + 3)
    semilogx(r_arr.*off,(ihl - ihlmid)./errmid,'.','color',...
        get_color(ifield-3),'markersize',10);hold on
    errorbar(r_arr.*off,(ihl - ihlmid)./errmid,ihl_err./errmid,'.','color',...
        get_color(ifield-3),'markersize',10);
    xlim([4e-1,1e3])
    xlabel('arcsec', 'fontsize',15)
    ylabel('scaled excess', 'fontsize',15)

end
end

savename=strcat(pltsavedir,'excess_all');
print(savename,'-dpng');%close
%%
figure
setwinsize(gcf,1200,300)

for im=1:3
ihldat = zeros([5,numel(ihl)]);
errdat = zeros([5,numel(ihl)]);
m_min = ihlprofdat.data(im).m_min;
m_max = ihlprofdat.data(im).m_max;

for ifield = 4:8
    dt=get_dark_times(flight,inst,ifield);

    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;

    ihl = excess(ifield).ihl(im).ihl;
    ihl_err = excess(ifield).ihl(im).err;
    ihldat(ifield-3,:) = ihl;
    errdat(ifield-3,:) = ihl_err;
end
ihltot = sum(ihldat./errdat.^2)./sum(1./errdat.^2);
ihltot_err = sqrt(1./(sum(1./errdat.^2)));

sp = find(r_arr > 100);
v = ihltot(sp);
e = ihltot_err(sp);
v = sum(v./e.^2)./sum(1./e.^2);
e = sqrt(1./sum(1./e.^2));

subplot(1,3,im)
loglog(r_arr.*off,ihltot,'k.','markersize',10);hold on
fill([r_arr(sp),flip(r_arr(sp))],[(v+e) * ones(size(sp)),(v-e) * ones(size(sp))],...
        [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
if v-e < 0
fill([r_arr(sp),flip(r_arr(sp))],[(v+e) * ones(size(sp)),1e-8 * ones(size(sp))],...
        [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');    
end
legend({'averaged exess','excess > 100 arcsec'})
legend boxoff

errorbar(r_arr.*off,ihltot,ihltot_err,'k.','markersize',10);
loglog(r_arr.*off,-ihltot,'ko','markersize',5);hold on
errorbar(r_arr.*off,-ihltot,ihltot_err,'ko','markersize',5);

ylim([1e-4,1e3])
xlim([4e-1,1e3])
title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);
xlabel('arcsec', 'fontsize',15);
ylabel('I [nW/m^2/sr]', 'fontsize',15);

end
savename=strcat(pltsavedir,'excess_avg');
print(savename,'-dpng');%close

return