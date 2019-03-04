function plot_ihl_excess(flight,inst)

mypaths=get_paths(flight);


pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

field_arr = 4:8;
for ifield=field_arr
    dt=get_dark_times(flight,inst,ifield);
    loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
    load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');
    excessdat.m_min_arr = ihlprofdat.m_min_arr;
    excessdat.m_max_arr = ihlprofdat.m_max_arr;
    excessdat.r_arr = ihlprofdat.r_arr;
    for im=10:12%1:numel(excessdat.m_min_arr)
        excessdat.field(ifield).diff(im).prof = ihlprofdat.excess(im).diff;
        excessdat.field(ifield).diff(im).err = ihlprofdat.excess(im).diff_err;
    end
end


%%% calculate the weighted average profile %%%
for im=10:12%1:numel(excessdat.m_min_arr)
    prof_tot = zeros(size(excessdat.r_arr));
    err_tot = zeros(size(excessdat.r_arr));
    for ifield=field_arr
        prof = excessdat.field(ifield).diff(im).prof;
        err = excessdat.field(ifield).diff(im).err;
        prof_tot = prof_tot + prof./err.^2;
        err_tot = err_tot + 1./err.^2;
    end
    prof_tot = prof_tot./err_tot;
    err_tot = sqrt(1./ err_tot);
    excessdat.avg(im).prof = prof_tot;
    excessdat.avg(im).err = err_tot;
end

savedir=strcat(mypaths.alldat,'TM',num2str(inst));
save(sprintf('%s/excessdat',savedir),'excessdat');

%%% plot the excess profile of each field %%%
figure
setwinsize(gcf,1200,300)

for im=10:12
    m_min = excessdat.m_min_arr(im);
    m_max = excessdat.m_max_arr(im);
    
    subplot(1,3,im-9)
    h=zeros(1,5);
    for ifield=field_arr
        dt=get_dark_times(flight,inst,ifield);
        r_arr = excessdat.r_arr;
        diff = excessdat.field(ifield).diff(im).prof;
        diff_err = excessdat.field(ifield).diff(im).err;
        
        h(ifield-3) = loglog(r_arr,diff,'.','color',get_color(ifield-3),...
            'markersize',10,'DisplayName',dt.name);hold on
        errorbar(r_arr,diff,diff_err,'.','color',get_color(ifield-3));
        loglog(r_arr,-diff,'o','color',get_color(ifield-3),...
            'markersize',5,'DisplayName',dt.name);
        errorbar(r_arr,-diff,diff_err,'o',...
            'color',get_color(ifield-3),'markersize',5);

    end
    xlim([4e-1,1e3])
    ylim([1e-4,1e0])
    xlabel('arcsec')
    ylabel('normalized <I_{stack}>')
    leg=legend(h,'Location','northeast');
    set(leg,'fontsize',10)
    legend boxoff
    title(strcat(num2str(m_min),'<mAB(y band)<',num2str(m_max)));

end
savename=strcat(pltsavedir,'rprof_excess_all');
print(savename,'-dpng');close

%%% plot the average profile %%%
figure
setwinsize(gcf,1200,300)

for im=10:12
    m_min = excessdat.m_min_arr(im);
    m_max = excessdat.m_max_arr(im);
    
    subplot(1,3,im-9)
    r_arr = excessdat.r_arr;
    diff = excessdat.avg(im).prof;
    diff_err = excessdat.avg(im).err;

    loglog(r_arr,diff,'.','markersize',10);hold on
    errorbar(r_arr,diff,diff_err,'k.');
    loglog(r_arr,-diff,'ko','markersize',5);
    errorbar(r_arr,-diff,diff_err,'ko','markersize',5);

    xlim([4e-1,1e3])
    ylim([1e-4,1e0])
    xlabel('arcsec')
    ylabel('normalized <I_{stack}>')
    title(strcat(num2str(m_min),'<mAB(y band)<',num2str(m_max)));

end
savename=strcat(pltsavedir,'rprof_excess_avg');
print(savename,'-dpng');close

return