flight = 40030;
inst = 1;
ifield = 8;
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/%s_mcerrdat_addnoise',loaddir,dt.name),'mcerrdat');
mcerrdatn = mcerrdat;
load(sprintf('%s/%s_mcerrdat',loaddir,dt.name),'mcerrdat');

Nsub = 10;

im = 3;    
figure
setwinsize(gcf,1200,300)
subplot(1,3,1)

r_arr = mcerrdat(im).dat.r_arr;
prof = mcerrdat(im).dat.profspserr_mat(1,:);
loglog(r_arr,prof,'r','DisplayName','\sigma_{annuli} / sqrt(N)');hold on
prof = mcerrdat(im).dat.profsps_std;
loglog(r_arr,prof,'b','linewidth', 3,'DisplayName','\sigma_{sims}');
loglog(r_arr,prof./sqrt(100),'b--','linewidth', 3,...
    'DisplayName','\sigma_{sims}/sqrt(100)');

xlim([4e-1,1e3])
ylim([1e-5,1e3])
h=legend('show','Location','southwest');
set(h,'fontsize',10)
legend boxoff
    
for isub = 1:Nsub
    prof = mcerrdat(im).dat.profspserr_mat(isub,:);
    loglog(r_arr,prof,'r');hold on
end

xlabel('arcsec', 'fontsize',20);
ylabel('\sigma_{stack}','fontsize',20);
title('Sim Map Stack');

subplot(1,3,2)

r_arr = mcerrdat(im).dat.r_arr;
prof = mcerrdat(im).dat.profscberr_mat(1,:);
loglog(r_arr,prof,'r','DisplayName','\sigma_{annuli} / sqrt(N)');hold on
prof = mcerrdat(im).dat.profscb_std;
loglog(r_arr,prof,'b','linewidth', 3,'DisplayName','\sigma_{sims}');
loglog(r_arr,prof./sqrt(100),'b--','linewidth', 3,...
    'DisplayName','\sigma_{sims}/sqrt(100)');

xlim([4e-1,1e3])
ylim([1e-5,1e3])
h=legend('show','Location','southwest');
set(h,'fontsize',10)
legend boxoff
    
for isub = 1:Nsub
    prof = mcerrdat(im).dat.profscberr_mat(isub,:);
    loglog(r_arr,prof,'r');hold on
end
xlabel('arcsec', 'fontsize',20);
ylabel('\sigma_{stack}','fontsize',20);
title('CIBER Map Stack');


subplot(1,3,3)

r_arr = mcerrdatn(im).dat.r_arr;
prof = mcerrdatn(im).dat.profspserr_mat(1,:);
loglog(r_arr,prof,'r','DisplayName','\sigma_{annuli} / sqrt(N)');hold on
prof = mcerrdatn(im).dat.profsps_std;
loglog(r_arr,prof,'b','linewidth', 3,'DisplayName','\sigma_{sims}');
loglog(r_arr,prof./sqrt(100),'b--','linewidth', 3,...
    'DisplayName','\sigma_{sims}/sqrt(100)');

xlim([4e-1,1e3])
ylim([1e-5,1e3])
h=legend('show','Location','southwest');
set(h,'fontsize',10)
legend boxoff
    
for isub = 1:Nsub
    prof = mcerrdatn(im).dat.profspserr_mat(isub,:);
    loglog(r_arr,prof,'r');hold on
end

xlabel('arcsec', 'fontsize',20);
ylabel('\sigma_{stack}','fontsize',20);
title('Sim Map + RN + PN Stack');


% subplot(1,3,2)
% prof = mcerrdatn(im).dat.profsps_std;
% loglog(r_arr,prof,'g','linewidth', 3,'DisplayName','\sigma_{sims}');


pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/DoF_test/'));
savename=strcat(pltsavedir,'data');
print(savename,'-dpng');%close