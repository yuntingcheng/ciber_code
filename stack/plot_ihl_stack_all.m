function plot_ihl_stack_all(flight,inst,ifield)

mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);


pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');

%%% write the basic info %%%
r_arr = ihlprofdat.r_arr;
psf_arr = ihlprofdat.psf_arr;
m_min = ihlprofdat.m_min_arr(10);
m_max = ihlprofdat.m_max_arr(12);
counts = sum(ihlprofdat.counts_arr(10:12));
countg = sum(ihlprofdat.countg_arr(10:12));

%%% plot the stacking profile %%%
    
figure
setwinsize(gcf,1000,800)

profscb = ihlprofdat.all.data.profscb;
profscb_err = ihlprofdat.all.data.profscb_err;
profsps = ihlprofdat.all.data.profsps;
profsps_err = ihlprofdat.all.data.profsps_err;
profgcb = ihlprofdat.all.data.profgcb;
profgcb_err = ihlprofdat.all.data.profgcb_err;
profgps = ihlprofdat.all.data.profgps;
profgps_err = ihlprofdat.all.data.profgps_err;

subplot(2,2,1)
errorbar(r_arr,profscb,profscb_err,'r.-','DisplayName','CBstars');hold on
errorbar(r_arr,profgcb,profgcb_err,'b.-','DisplayName','CBgals');
errorbar(r_arr,profsps,profsps_err,'m.-','DisplayName','PSstars');
errorbar(r_arr,profgps,profgps_err,'c.-','DisplayName','PSgals');     

xlabel('arcsec')
ylabel('<I_{stack}>')
title('stack unnormalized (linear)',...
    'interpreter','latex','fontsize',15)
xlim([4e-1,7e2])
set(gca, 'XScale', 'log')
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff
    
subplot(2,2,2)
loglog(r_arr.*0.98,profscb,'r.','markersize',10);hold on
loglog(r_arr.*1.02,profgcb,'b.','markersize',10);
errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
loglog(r_arr.*0.98,-profscb,'ro','markersize',5);hold on
loglog(r_arr.*1.02,-profgcb,'bo','markersize',5);
errorbar(r_arr.*0.98,-profscb,profscb_err,'ro','markersize',5);
errorbar(r_arr.*1.02,-profgcb,profgcb_err,'bo','markersize',5);

loglog(r_arr.*0.98,profsps,'m.','markersize',10);hold on
loglog(r_arr.*1.02,profgps,'c.','markersize',10);
errorbar(r_arr.*0.98,profsps,profsps_err,'m.','markersize',10);
errorbar(r_arr.*1.02,profgps,profgps_err,'c.','markersize',10);
loglog(r_arr.*0.98,-profsps,'mo','markersize',5);hold on
loglog(r_arr.*1.02,-profgps,'co','markersize',5);
errorbar(r_arr.*0.98,-profsps,profsps_err,'mo','markersize',5);
errorbar(r_arr.*1.02,-profgps,profgps_err,'co','markersize',5);

xlim([4e-1,1e3])
xlabel('arcsec')
ylabel('<I_{stack}>')
title('stack unnormalized (log)',...
    'interpreter','latex','fontsize',15)

profscb = ihlprofdat.all.norm.profscb;
profscb_err = ihlprofdat.all.norm.profscb_err;
profsps = ihlprofdat.all.norm.profsps;
profsps_err = ihlprofdat.all.norm.profsps_err;
profgcb = ihlprofdat.all.norm.profgcb;
profgcb_err = ihlprofdat.all.norm.profgcb_err;
profgps = ihlprofdat.all.norm.profgps;
profgps_err = ihlprofdat.all.norm.profgps_err;

subplot(2,2,3)
errorbar(r_arr,profscb,profscb_err,'r.-','DisplayName','CBstars');hold on
errorbar(r_arr,profgcb,profgcb_err,'b.-','DisplayName','CBgals');
errorbar(r_arr,profsps,profsps_err,'m.-','DisplayName','PSstars');
errorbar(r_arr,profgps,profgps_err,'c.-','DisplayName','PSgals');     

xlabel('arcsec')
ylabel('<I_{stack}>')
title('stack unnormalized (linear)',...
    'interpreter','latex','fontsize',15)
xlim([4e-1,7e2])
set(gca, 'XScale', 'log')
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff

subplot(2,2,4)
loglog(r_arr.*0.98,profscb,'r.','markersize',10);hold on
loglog(r_arr.*1.02,profgcb,'b.','markersize',10);
errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
loglog(r_arr.*0.98,-profscb,'ro','markersize',5);hold on
loglog(r_arr.*1.02,-profgcb,'bo','markersize',5);
errorbar(r_arr.*0.98,-profscb,profscb_err,'ro','markersize',5);
errorbar(r_arr.*1.02,-profgcb,profgcb_err,'bo','markersize',5);

loglog(r_arr.*0.98,profsps,'m.','markersize',10);hold on
loglog(r_arr.*1.02,profgps,'c.','markersize',10);
errorbar(r_arr.*0.98,profsps,profsps_err,'m.','markersize',10);
errorbar(r_arr.*1.02,profgps,profgps_err,'c.','markersize',10);
loglog(r_arr.*0.98,-profsps,'mo','markersize',5);hold on
loglog(r_arr.*1.02,-profgps,'co','markersize',5);
errorbar(r_arr.*0.98,-profsps,profsps_err,'mo','markersize',5);
errorbar(r_arr.*1.02,-profgps,profgps_err,'co','markersize',5);


loglog(r_arr,psf_arr,'k--');
xlim([4e-1,1e3])
ymin = [abs(profscb),abs(profgcb),abs(profsps),abs(profgps)];
ymin=ymin(find(ymin>0));
ymin=floor(log10(min(ymin)));
ymin=10^ymin;
ylim([ymin,2e0])
xlabel('arcsec')
ylabel('<I_{stack}>')
title('stack unnormalized (log)',...
    'interpreter','latex','fontsize',15)

suptitle(strcat('PanSTARRS',{' '},...
    num2str(m_min),'<mAB(y band)<',num2str(m_max),...
    '(',num2str(counts),' stars, ',num2str(countg), ' gals)'));

savename=strcat(pltsavedir,...
    dt.name,'_rprof',num2str(m_min),'_',num2str(m_max),'bksub');
print(savename,'-dpng');%close

%%% plot the excess profile %%%
figure

diffcb = ihlprofdat.all.excess.diffcb;
diffcb_err = ihlprofdat.all.excess.diffcb_err;
diffps = ihlprofdat.all.excess.diffps;
diffps_err = ihlprofdat.all.excess.diffsp_err;
diff = ihlprofdat.all.excess.diff;
diff_err = ihlprofdat.all.excess.diff_err;

loglog(r_arr.*0.98,diffcb,'r.','markersize',10,...
    'DisplayName','CB diff = CB gal - CB stars');hold on
loglog(r_arr.*1.02,diffps,'b.','markersize',10,...
    'DisplayName','PS diff = PS gal - PS stars');
loglog(r_arr,diff,'k.','markersize',10,...
    'DisplayName','CB diff - PS diff');
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff    
errorbar(r_arr.*0.98,diffcb,diffcb_err,'r.');
errorbar(r_arr.*1.02,diffps,diffps_err,'b.');
errorbar(r_arr,diff,diff_err,'k.');


loglog(r_arr.*0.98,-diffcb,'ro','markersize',5);hold on
loglog(r_arr.*1.02,-diffps,'bo','markersize',5);
loglog(r_arr,-diff,'ko','markersize',5);
errorbar(r_arr.*0.98,-diffcb,diffcb_err,'ro','markersize',5);
errorbar(r_arr.*1.02,-diffps,diffps_err,'bo','markersize',5);
errorbar(r_arr,-diff,diff_err,'ko','markersize',5);

xlim([4e-1,1e3])
xlabel('arcsec')
ylabel('normalized <I_{stack}>')

title(strcat(num2str(m_min),'<mAB(y band)<',num2str(m_max)));

savename=strcat(pltsavedir,dt.name,'_rprof_excess_all');
print(savename,'-dpng');%close

return