function plot_ihl_stack_old1(flight,inst)
%%
mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');

%%% plot the stacking profile %%%
figure
setwinsize(gcf,1200,300)
for im=1:3
    m_min = ihlprofdat(im).m_min;
    m_max = ihlprofdat(im).m_max;
    counts = ihlprofdat(im).counts;
    countg = ihlprofdat(im).countg;
    r_arr = ihlprofdat.r_arr;

    profscb = ihlprofdat.data(im).profscb;
    profscb_err = ihlprofdat.data(im).profscb_err;
    profsps = ihlprofdat.data(im).profsps;
    profsps_err = ihlprofdat.data(im).profsps_err;
    profgcb = ihlprofdat.data(im).profgcb;
    profgcb_err = ihlprofdat.data(im).profgcb_err;
    profgps = ihlprofdat.data(im).profgps;
    profgps_err = ihlprofdat.data(im).profgps_err;
    
    subplot(1,3,im)
    loglog(r_arr.*0.98,profscb,'r.','markersize',10);hold on
    loglog(r_arr.*1.02,profgcb,'b.','markersize',10);
    loglog(r_arr.*0.98,profsps,'m.','markersize',10);hold on
    loglog(r_arr.*1.02,profgps,'c.','markersize',10);
    h=legend({'CIBER stars','CIBER galaxies','Sim stars','Sim galaxies'},...
        'Location','northeast');
    set(h,'fontsize',10)
    legend boxoff
    
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    loglog(r_arr.*0.98,-profscb,'ro','markersize',5);hold on
    loglog(r_arr.*1.02,-profgcb,'bo','markersize',5);
    errorbar(r_arr.*0.98,-profscb,profscb_err,'ro','markersize',5);
    errorbar(r_arr.*1.02,-profgcb,profgcb_err,'bo','markersize',5);
    errorbar(r_arr.*0.98,profsps,profsps_err,'m.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'c.','markersize',10);
    loglog(r_arr.*0.98,-profsps,'mo','markersize',5);hold on
    loglog(r_arr.*1.02,-profgps,'co','markersize',5);
    errorbar(r_arr.*0.98,-profsps,profsps_err,'mo','markersize',5);
    errorbar(r_arr.*1.02,-profgps,profgps_err,'co','markersize',5);

    xlim([4e-1,1e3])
    ylim([1e-2,1e4])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
      ' (',num2str(counts),' stars, ',num2str(countg), ' galaxies)'),'fontsize',15)
end
savename=strcat(pltsavedir,dt.name,'_rprof');
print(savename,'-dpng');close
%%
figure
setwinsize(gcf,1200,300)
for im=1:3
    m_min = ihlprofdat(im).m_min;
    m_max = ihlprofdat(im).m_max;
    counts = ihlprofdat(im).counts;
    countg = ihlprofdat(im).countg;
    r_arr = ihlprofdat.r_arr;
    psf_arr = ihlprofdat.psf_arr;

    profscb = ihlprofdat.norm(im).profscb;
    profscb_err = ihlprofdat.norm(im).profscb_err;
    profsps = ihlprofdat.norm(im).profsps;
    profsps_err = ihlprofdat.norm(im).profsps_err;
    profgcb = ihlprofdat.norm(im).profgcb;
    profgcb_err = ihlprofdat.norm(im).profgcb_err;
    profgps = ihlprofdat.norm(im).profgps;
    profgps_err = ihlprofdat.norm(im).profgps_err;

    subplot(1,3,im)
    loglog(r_arr.*0.98,profscb,'r.','markersize',10);hold on
    loglog(r_arr.*1.02,profgcb,'b.','markersize',10);
    loglog(r_arr.*0.98,profsps,'m.','markersize',10);hold on
    loglog(r_arr.*1.02,profgps,'c.','markersize',10);
    loglog(r_arr, psf_arr,'k--')
    h=legend({'CIBER stars','CIBER galaxies','Sim stars','Sim galaxies','PSF'},...
        'Location','southwest');
    set(h,'fontsize',10)
    legend boxoff
    
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    loglog(r_arr.*0.98,-profscb,'ro','markersize',5);hold on
    loglog(r_arr.*1.02,-profgcb,'bo','markersize',5);
    errorbar(r_arr.*0.98,-profscb,profscb_err,'ro','markersize',5);
    errorbar(r_arr.*1.02,-profgcb,profgcb_err,'bo','markersize',5);
    errorbar(r_arr.*0.98,profsps,profsps_err,'m.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'c.','markersize',10);
    loglog(r_arr.*0.98,-profsps,'mo','markersize',5);hold on
    loglog(r_arr.*1.02,-profgps,'co','markersize',5);
    errorbar(r_arr.*0.98,-profsps,profsps_err,'mo','markersize',5);
    errorbar(r_arr.*1.02,-profgps,profgps_err,'co','markersize',5);

    xlim([4e-1,1e3])
    ylim([1e-6,1.2])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I / I(r=0)', 'fontsize',15)
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
      ' (',num2str(counts),' stars, ',num2str(countg), ' galaxies)'),'fontsize',15)
end
savename=strcat(pltsavedir,dt.name,'_rprof_norm');
print(savename,'-dpng');close
%%
figure
setwinsize(gcf,1200,300)

for im=1:3
    m_min = ihlprofdat(im).m_min;
    m_max = ihlprofdat(im).m_max;
    counts = ihlprofdat(im).counts;
    countg = ihlprofdat(im).countg;
    r_arr = ihlprofdat.r_arr;
    
    subplot(1,3,im)
    diff = ihlprofdat.excess(im).diff;
    loglog(r_arr,diff,'k.','markersize',10);hold on
    diff = ihlprofdat.excess(im).diffps;
    loglog(r_arr,diff,'b.','markersize',10);hold on
    diff = ihlprofdat.excess(im).diffcb;
    loglog(r_arr,diff,'r.','markersize',10);hold on
    h=legend({'CIBER gals - CIBER stars','Sim gals - Sim stars'...
        'CIBER excess - Sim excess'},...
        'Location','northeast');
    set(h,'fontsize',10)
    legend boxoff
  
    diff = ihlprofdat.excess(im).diff;
    diff_err = ihlprofdat.excess(im).diff_err;
    errorbar(r_arr,diff,diff_err,'k.');
    loglog(r_arr,-diff,'ko','markersize',5);
    errorbar(r_arr,-diff,diff_err,'ko','markersize',5);
    diff = ihlprofdat.excess(im).diffcb;
    diff_err = ihlprofdat.excess(im).diffcb_err;
    errorbar(r_arr,diff,diff_err,'r.');
    loglog(r_arr,-diff,'ro','markersize',5);
    errorbar(r_arr,-diff,diff_err,'ro','markersize',5);
    diff = ihlprofdat.excess(im).diffps;
    diff_err = ihlprofdat.excess(im).diffps_err;
    errorbar(r_arr,diff,diff_err,'b.');
    loglog(r_arr,-diff,'bo','markersize',5);
    errorbar(r_arr,-diff,diff_err,'bo','markersize',5);

    xlim([4e-1,1e3])
    ylim([1e-5,1e0])
    xlabel('arcsec','fontsize',15)
    ylabel('I / I(r=0) excess','fontsize',15)
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
      ' (',num2str(counts),' stars, ',num2str(countg), ' galaxies)'),'fontsize',15)
end
suptitle(dt.name);

savename=strcat(pltsavedir,dt.name,'_rprof_excess');
print(savename,'-dpng');%close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot excess of all fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
setwinsize(gcf,1200,600)
xplot_off = [0,0,0,0.97,0.98,1,1.02,1.03];
mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

for im=1:3
for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');
%%% write the basic info %%%
xoff = xplot_off(ifield);
m_min = ihlprofdat(im).m_min;
m_max = ihlprofdat(im).m_max;
r_arr = ihlprofdat.r_arr;

diff = ihlprofdat.excess(im).diff;

subplot(2,3,im)
loglog(r_arr.*xoff,diff,'.','color',get_color(ifield-3),...
    'markersize',10,'DisplayName',dt.name);hold on
subplot(2,3,im+3)
semilogx(r_arr.*xoff,diff,'.','color',get_color(ifield-3),...
    'markersize',10,'DisplayName',dt.name);hold on
end
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff
end


for im=1:3
for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');
xoff = xplot_off(ifield);
m_min = ihlprofdat(im).m_min;
m_max = ihlprofdat(im).m_max;
r_arr = ihlprofdat.r_arr;

diff = ihlprofdat.excess(im).diff;
diff_err = ihlprofdat.excess(im).diff_err;

subplot(2,3,im)
errorbar(r_arr.*xoff,diff,diff_err,'.','color',get_color(ifield-3));
loglog(r_arr.*xoff,-diff,'o','color',get_color(ifield-3),'markersize',5);
errorbar(r_arr.*xoff,-diff,diff_err,'o','color',get_color(ifield-3),...
    'markersize',5);
subplot(2,3,im+3)
errorbar(r_arr.*xoff,diff,diff_err,'.','color',get_color(ifield-3));
end

subplot(2,3,im)
xlim([4e-1,1e3])
ylim([1e-5,1e0])
xlabel('arcsec','fontsize',15)
ylabel('I / I(r=0) excess','fontsize',15)
title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15)
subplot(2,3,im+3)
xlim([4e-1,1e3])
xlabel('arcsec','fontsize',15)
ylabel('I / I(r=0) excess','fontsize',15)
title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15)
end
savename=strcat(pltsavedir,'rprof_excess_all');
print(savename,'-dpng');%close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot excess of weighted avg of all fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
setwinsize(gcf,1200,600)
mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

for im=1:3

    m_min = ihlprofdat(im).m_min;
    m_max = ihlprofdat(im).m_max;
    mu = 0;
    sig = 0;
    for ifield = 4:8
        dt=get_dark_times(flight,inst,ifield);
        loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
        load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');
        diff = ihlprofdat.excess(im).diff;
        diff_err = ihlprofdat.excess(im).diff_err;
        mu = mu + diff./diff_err.^2;
        sig = sig + 1./diff_err.^2;
    end
    mu = mu ./ sig;
    sig = sqrt(1./sig);
    
    subplot(2,3,im)
    loglog(r_arr,mu,'k.','markersize',10);hold on
    errorbar(r_arr,mu,sig,'k.');
    loglog(r_arr,-mu,'ko','markersize',5);
    errorbar(r_arr,-mu,sig,'ko','markersize',5);
    xlim([4e-1,1e3])
    ylim([1e-5,1e0])
    xlabel('arcsec','fontsize',15)
    ylabel('I / I(r=0) excess','fontsize',15)
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15)

    subplot(2,3,im+3)
    semilogx(r_arr,mu,'k.','markersize',10);hold on
    errorbar(r_arr,mu,sig,'k.');
    xlim([4e-1,1e3])
    xlabel('arcsec','fontsize',15)
    ylabel('I / I(r=0) excess','fontsize',15)
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15)

end
savename=strcat(pltsavedir,'rprof_excess_avg');
print(savename,'-dpng');%close

return