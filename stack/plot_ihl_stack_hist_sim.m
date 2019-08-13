function plot_ihl_stack_hist_sim(flight,inst,ifield,rvir)
%%
mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));
dt=get_dark_times(flight,inst,ifield);
%%
for f_ihl=[0,0.1,0.5]

%%% plot the stacking profile %%%
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if rvir==1 || f_ihl==0
    load(sprintf('%s/%s_ihlprofdatsim_hist%d',...
        loaddir,dt.name,f_ihl*100),'ihlprofdat');
else
    load(sprintf('%s/%s_ihlprofdatsim_hist%d_rv%d',...
        loaddir,dt.name,f_ihl*100,rvir),'ihlprofdat');
end

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

    profscbbk = ihlprofdat.bk(im).profscb;
    profscb_errbk = ihlprofdat.bk(im).profscb_err;
    profspsbk = ihlprofdat.bk(im).profsps;
    profsps_errbk = ihlprofdat.bk(im).profsps_err;

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
    plot(r_arr,profscbbk,'k','markersize',10);
    h=legend({'CIBER stars','CIBER galaxies','background'},...
        'Location','northwest');
    set(h,'fontsize',10)
    legend boxoff
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    errorbar(r_arr,profscbbk,profscb_errbk,'k.','markersize',10);
    xlim([4e-1,1e3])
    ylim([-1,10])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)

    subplot(3,3,im + 6)
    semilogx(r_arr.*0.98,profsps,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgps,'b','markersize',10);
    plot(r_arr,profspsbk,'k','markersize',10);
    h=legend({'Sim stars','Sim galaxies'},...
        'Location','northwest');
    set(h,'fontsize',10)
    legend boxoff
    errorbar(r_arr.*0.98,profsps,profsps_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'b.','markersize',10);
    errorbar(r_arr,profspsbk,profsps_errbk,'k.','markersize',10);
    xlim([4e-1,1e3])
    ylim([-0.05,0.1])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end
s=suptitle(strcat(dt.name,'  f_{IHL}=',num2str(f_ihl)));
set(s,'FontSize',15);
savename=strcat(pltsavedir,dt.name,'_stackprof_hist_sim',num2str(f_ihl*100));
% print(savename,'-dpng');%close

%%% plot the excess profile%%%

dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if rvir==1 || f_ihl==0
    load(sprintf('%s/%s_ihlprofdatsim_hist%d',...
        loaddir,dt.name,f_ihl*100),'ihlprofdat');
else
    load(sprintf('%s/%s_ihlprofdatsim_hist%d_rv%d',...
        loaddir,dt.name,f_ihl*100,rvir),'ihlprofdat');
end

figure
setwinsize(gcf,1200,360)
for im=1:3
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
    
    subplot(1,3,im)
    loglog(r_arr.*0.98,ecb,'r.','markersize',10);hold on
    loglog(r_arr.*1.02,eps,'b.','markersize',10);hold on
    loglog(r_arr,ihl,'k.','markersize',10);hold on
    
    sp = find(r_arr > 100);
    v = ihl(sp);
    e = ihl_err(sp);
    v = sum(v./e.^2)./sum(1./e.^2);
    e = sqrt(1./sum(1./e.^2));
    
    if v-e >=0
    fill([r_arr(sp),flip(r_arr(sp))],...
        [(v+e) * ones(size(sp)),(v-e) * ones(size(sp))],...
            [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
    else
    fill([r_arr(sp),flip(r_arr(sp))],...
        [(v+e) * ones(size(sp)),1e-8 * ones(size(sp))],...
            [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');    
    end

    h=legend({'CIBER gals - stars','Sim gals - stars', ...
        'CIBER excess - Sim excess','excess > 100 arcsec'},'Location','northeast');
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
s=suptitle(strcat(dt.name,'  f_{IHL}=',num2str(f_ihl)));
set(s,'FontSize',15);
savename=strcat(pltsavedir,dt.name,'_excessprof_hist_sim',num2str(f_ihl*100));
print(savename,'-dpng');%close
%%%
figure

for im=1:3
    off = 1 + (im - 2).*0.02;
    r_arr = ihlprofdat.r_arr;
    ihl = ihlprofdat.excess(im).diff;
    loglog(r_arr.*off,ihl,'.','color',get_color(im),'markersize',10);hold on
end
legend({'16<m<17','17<m<18', '18<m<19'});
legend boxoff

for im=1:3
    off = 1 + (im - 2).*0.02;
    r_arr = ihlprofdat.r_arr;
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;
    errorbar(r_arr.*off,ihl,ihl_err,'.','color',get_color(im),'markersize',10);
    loglog(r_arr.*off,-ihl,'o','color',get_color(im),'markersize',5);hold on
    errorbar(r_arr.*off,-ihl,ihl_err,'o','color',get_color(im),'markersize',5);
end

ylim([1e-4,1e3])
xlim([4e-1,1e3])
xlabel('arcsec', 'fontsize',15);
ylabel('I [nW/m^2/sr]', 'fontsize',15);

s=title(strcat(dt.name,'  f_{IHL}=',num2str(f_ihl)));
set(s,'FontSize',15);

savename=strcat(pltsavedir,'excess_avg_all_hist_sim',num2str(f_ihl*100));
% print(savename,'-dpng');%close
end
%%
figure
setwinsize(gcf,1200,360)
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
c_arr = {'k','b','r','m'};
off_arr=[1,0.98,0.96,1.02];
count=0;
for f_ihl=[-999,0,0.1,0.5]
count=count+1;
if f_ihl==-999
    load(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');
else
    if rvir==1 || f_ihl==0
        load(sprintf('%s/%s_ihlprofdatsim_hist%d',...
            loaddir,dt.name,f_ihl*100),'ihlprofdat');
    else
        load(sprintf('%s/%s_ihlprofdatsim_hist%d_rv%d',...
            loaddir,dt.name,f_ihl*100,rvir),'ihlprofdat');
    end

end
c = c_arr{count};
off = off_arr(count);
for im=1:3
    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;
    r_arr = ihlprofdat.r_arr;
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;
    
    subplot(1,3,im)
    hh(count)=loglog(r_arr*off,ihl,'.','color',c,'markersize',10);hold on
    

    errorbar(r_arr*off,ihl,ihl_err,'.','color',c,'markersize',10);
    loglog(r_arr*off,-ihl,'o','color',c,'markersize',5);hold on
    errorbar(r_arr*off,-ihl,ihl_err,'o','color',c,'markersize',5);
    
    xlim([4e-1,1e3])
    ylim([1e-4,1e4])
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15)
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end
end

for im=1:3
    subplot(1,3,im)
    h=legend(hh,'CIBER data','Sim no IHL','Sim 10 % IHL','Sim 50 % IHL');
    set(h,'fontsize',10,'Location','northeast')
    legend boxoff
end

savename=strcat(pltsavedir,dt.name,'_excessprof_hist_sim');
% print(savename,'-dpng');%close
%%
figure
setwinsize(gcf,1200,600)
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
c_arr = {'k','b','r','c','m'};
off_arr=[1,0.98,0.96,1.02,1.04];

f_ihl=0.1;
count=0;
for rvir=[-999,0,1,2,3]
count=count+1;
if rvir==-999
    load(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');
else
    if rvir==0
        load(sprintf('%s/%s_ihlprofdatsim_hist%d',...
            loaddir,dt.name,0),'ihlprofdat');
    elseif rvir==1
        load(sprintf('%s/%s_ihlprofdatsim_hist%d',...
            loaddir,dt.name,f_ihl*100),'ihlprofdat');
    else
        load(sprintf('%s/%s_ihlprofdatsim_hist%d_rv%d',...
            loaddir,dt.name,f_ihl*100,rvir),'ihlprofdat');
    end

end
c = c_arr{count};
off = off_arr(count);
for im=1:3
    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;
    r_arr = ihlprofdat.r_arr;
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;
    
    subplot(2,3,im)
    hh(count)=loglog(r_arr*off,ihl,'.','color',c,'markersize',10);hold on
    

    errorbar(r_arr*off,ihl,ihl_err,'.','color',c,'markersize',10);
    loglog(r_arr*off,-ihl,'o','color',c,'markersize',5);hold on
    errorbar(r_arr*off,-ihl,ihl_err,'o','color',c,'markersize',5);
    
    xlim([4e-1,1e3])
    ylim([1e-4,1e4])
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15)
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end
end

for im=1:3
    subplot(2,3,im)
    h=legend(hh,'CIBER data','Sim no IHL',...
        '10 % IHL 1Rvir','10 % IHL 2Rvir','10 % IHL 3Rvir');
    set(h,'fontsize',10,'Location','northeast')
    legend boxoff
end

f_ihl=0.5;
count=0;
for rvir=[-999,0,1,2,3]
count=count+1;
if rvir==-999
    load(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');
else
    if rvir==0
        load(sprintf('%s/%s_ihlprofdatsim_hist%d',...
            loaddir,dt.name,0),'ihlprofdat');
    elseif rvir==1
        load(sprintf('%s/%s_ihlprofdatsim_hist%d',...
            loaddir,dt.name,f_ihl*100),'ihlprofdat');
    else
        load(sprintf('%s/%s_ihlprofdatsim_hist%d_rv%d',...
            loaddir,dt.name,f_ihl*100,rvir),'ihlprofdat');
    end

end
c = c_arr{count};
off = off_arr(count);
for im=1:3
    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;
    r_arr = ihlprofdat.r_arr;
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;
    
    subplot(2,3,im+3)
    hh(count)=loglog(r_arr*off,ihl,'.','color',c,'markersize',10);hold on
    

    errorbar(r_arr*off,ihl,ihl_err,'.','color',c,'markersize',10);
    loglog(r_arr*off,-ihl,'o','color',c,'markersize',5);hold on
    errorbar(r_arr*off,-ihl,ihl_err,'o','color',c,'markersize',5);
    
    xlim([4e-1,1e3])
    ylim([1e-4,1e4])
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15)
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end
end

for im=1:3
    subplot(2,3,im+3)
    h=legend(hh,'CIBER data','Sim no IHL',...
        '50 % IHL 1Rvir','50 % IHL 2Rvir','50 % IHL 3Rvir');
    set(h,'fontsize',10,'Location','northeast')
    legend boxoff
end



savename=strcat(pltsavedir,dt.name,'_excessprof_hist_sim_allrv');
print(savename,'-dpng');%close
%%
return