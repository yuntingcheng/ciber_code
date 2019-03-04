function plot_ihl_excess_all(flight,inst)
%%
mypaths=get_paths(flight);


pltsavedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

savedir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/excessdat',savedir),'excessdat');

m_min = excessdat.m_min_arr(10);
m_max = excessdat.m_max_arr(12);

for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
    load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');
    excessdat.all_field(ifield).diff.prof = ihlprofdat.all.excess.diff;
    excessdat.all_field(ifield).diff.err = ihlprofdat.all.excess.diff_err;
end


%%% calculate the weighted average profile %%%
prof_tot = zeros(size(excessdat.r_arr));
err_tot = zeros(size(excessdat.r_arr));
for ifield=4:8
    prof = excessdat.all_field(ifield).diff.prof;
    err = excessdat.all_field(ifield).diff.err;
    prof_tot = prof_tot + prof./err.^2;
    err_tot = err_tot + 1./err.^2;
end
prof_tot = prof_tot./err_tot;
err_tot = sqrt(1./ err_tot);
excessdat.all_avg.prof = prof_tot;
excessdat.all_avg.err = err_tot;

prof_tot = zeros(size(excessdat.r_arr));
err_tot = zeros(size(excessdat.r_arr));
for ifield=7:8
    prof = excessdat.all_field(ifield).diff.prof;
    err = excessdat.all_field(ifield).diff.err;
    prof_tot = prof_tot + prof./err.^2;
    err_tot = err_tot + 1./err.^2;
end
prof_tot = prof_tot./err_tot;
err_tot = sqrt(1./ err_tot);
excessdat.all_avg.prof78 = prof_tot;
excessdat.all_avg.err78 = err_tot;

savedir=strcat(mypaths.alldat,'TM',num2str(inst));
save(sprintf('%s/excessdat',savedir),'excessdat');

%%% plot the excess profile of each field %%%

figure
h=zeros(1,numel(4:8));
fieldcount = 0;
xplot_off = [0,0,0,0.97,0.98,1,1.02,1.03];
for ifield=4:8
    fieldcount = fieldcount + 1;
    dt = get_dark_times(flight,inst,ifield);
    r_arr = excessdat.r_arr;
    diff = excessdat.all_field(ifield).diff.prof;
    diff_err = excessdat.all_field(ifield).diff.err;

    xoff = xplot_off(ifield);
    h(fieldcount) = loglog(r_arr.*xoff,diff,'.','color',get_color(ifield-3),...
        'markersize',10,'DisplayName',dt.name);hold on
    errorbar(r_arr.*xoff,diff,diff_err,'.','color',get_color(ifield-3));
    loglog(r_arr.*xoff,-diff,'o','color',get_color(ifield-3),...
        'markersize',5,'DisplayName',dt.name);
    errorbar(r_arr.*xoff,-diff,diff_err,'o',...
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

savename=strcat(pltsavedir,'rprof_excess_allmag');
print(savename,'-dpng');

%%% plot the average profile %%%
figure
r_arr = excessdat.r_arr;
diff = excessdat.all_avg.prof;
diff_err = excessdat.all_avg.err;

loglog(r_arr,diff,'.','markersize',10);hold on
errorbar(r_arr,diff,diff_err,'k.');
loglog(r_arr,-diff,'ko','markersize',5);
errorbar(r_arr,-diff,diff_err,'ko','markersize',5);

xlim([4e-1,1e3])
ylim([1e-4,1e-1])
xlabel('arcsec')
ylabel('normalized <I_{stack}>')
title(strcat(num2str(m_min),'<mAB(y band)<',num2str(m_max)));

savename=strcat(pltsavedir,'rprof_excess_avg_all');
print(savename,'-dpng');

%%% plot only avg of field7 & 8 %%%
% figure
% r_arr = excessdat.r_arr;
% diff = excessdat.all_avg.prof78;
% diff_err = excessdat.all_avg.err78;
% 
% loglog(r_arr,diff,'.','markersize',10);hold on
% errorbar(r_arr,diff,diff_err,'k.');
% loglog(r_arr,-diff,'ko','markersize',5);
% errorbar(r_arr,-diff,diff_err,'ko','markersize',5);
% 
% xlim([4e-1,1e3])
% ylim([1e-4,1e-1])
% xlabel('arcsec')
% ylabel('normalized <I_{stack}>')
% title(strcat(num2str(m_min),'<mAB(y band)<',num2str(m_max)));
% 
% savename=strcat(pltsavedir,'rprof_excess_avg_all78');
% print(savename,'-dpng');

return