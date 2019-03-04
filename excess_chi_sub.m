flight = 40030;
inst = 1;
ifield = 8;

pltsavedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

stackdir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/%s_mcerrdat',loaddir,dt.name),'mcerrdat');

nsim = 100;

m_min_arr = stackmapdat(ifield).m_min_arr(10:12);
m_max_arr = stackmapdat(ifield).m_max_arr(10:12);
counts_arr = stackmapdat(ifield).count_stacks_arr(10:12);
countg_arr = stackmapdat(ifield).count_stackg_arr(10:12);

%%% write the basic info %%%
ihlprofdat.m_min_arr = m_min_arr;
ihlprofdat.m_max_arr = m_max_arr;
ihlprofdat.counts_arr = counts_arr;
ihlprofdat.countg_arr = countg_arr;

%%
im = 3;
dat = mcerrdat(im).dat;
r_arr = dat.r_arr;

m_min = m_min_arr(im);
m_max = m_max_arr(im);

rmask = 0;
for i = 4:8
    rmask = rmask + get_mask_radius(inst,i,(m_min + m_max)/2);
end
rmask = rmask / 5;
sp1 = find(r_arr > rmask & r_arr <= 100);
sp2 = find(r_arr > 100);


figure
setwinsize(gcf,1200,400)

subplot(2,4,1)
mu = mean(dat.profscb_mat);
sig = nanstd(dat.profscb_mat);
chi0_arr = zeros([1, dat.Nsubs]);
chi1_arr = zeros([1, dat.Nsubs]);
chi2_arr = zeros([1, dat.Nsubs]);
for i = 1:dat.Nsubs
    semilogx(r_arr, dat.profscb_mat(i,:)); hold on
    chi0_arr(i) = sum((dat.profscb_mat(i,:) - mu).^2 ./ sig.^2);
    chi1_arr(i) = sum((dat.profscb_mat(i,sp1) - mu(sp1)).^2 ./ sig(sp1).^2);
    chi2_arr(i) = sum((dat.profscb_mat(i,sp2) - mu(sp2)).^2 ./ sig(sp2).^2);
end
semilogx(r_arr, mu,'k--', 'linewidth',3); hold on
xlim([4e-1,1e3])
xlabel('arcsec', 'fontsize',15);
ylabel('<I>_{stack}','fontsize',15);
title('CIBER stars','fontsize', 15);
%legend('show');
subplot(2,4,5)
plot(chi0_arr./numel(r_arr),'k.','markersize',20);hold on
plot(chi1_arr./numel(sp1),'r.','markersize',20);
plot(chi2_arr./numel(sp2),'b.','markersize',20);
xlim([0.8,10.2])
ylabel('\chi^2 / N', 'fontsize', 15)
xlabel('sub-stack', 'fontsize', 15)

subplot(2,4,2)
mu = mean(dat.profsps_mat);
sig = nanstd(dat.profsps_mat);
chi0_arr = zeros([1, dat.Nsubs]);
chi1_arr = zeros([1, dat.Nsubs]);
chi2_arr = zeros([1, dat.Nsubs]);
for i = 1:dat.Nsubs
    semilogx(r_arr, dat.profsps_mat(i,:)); hold on
    chi0_arr(i) = sum((dat.profsps_mat(i,:) - mu).^2 ./ sig.^2);
    chi1_arr(i) = sum((dat.profsps_mat(i,sp1) - mu(sp1)).^2 ./ sig(sp1).^2);
    chi2_arr(i) = sum((dat.profsps_mat(i,sp2) - mu(sp2)).^2 ./ sig(sp2).^2);
end
semilogx(r_arr, mu,'k--', 'linewidth',3); hold on
xlim([4e-1,1e3])
xlabel('arcsec', 'fontsize',15);
ylabel('<I>_{stack}','fontsize',15);
title('Sim stars','fontsize', 15)
subplot(2,4,6)
plot(chi0_arr./numel(r_arr),'k.','markersize',20);hold on
plot(chi1_arr./numel(sp1),'r.','markersize',20);
plot(chi2_arr./numel(sp2),'b.','markersize',20);
xlim([0.8,10.2])
ylabel('\chi^2 / N', 'fontsize', 15)
xlabel('sub-stack', 'fontsize', 15)

subplot(2,4,3)
mu = mean(dat.profgcb_mat);
sig = nanstd(dat.profgcb_mat);
chi0_arr = zeros([1, dat.Nsubs]);
chi1_arr = zeros([1, dat.Nsubs]);
chi2_arr = zeros([1, dat.Nsubs]);
for i = 1:dat.Nsubs
    semilogx(r_arr, dat.profgcb_mat(i,:)); hold on
    chi0_arr(i) = sum((dat.profgcb_mat(i,:) - mu).^2 ./ sig.^2);
    chi1_arr(i) = sum((dat.profgcb_mat(i,sp1) - mu(sp1)).^2 ./ sig(sp1).^2);
    chi2_arr(i) = sum((dat.profgcb_mat(i,sp2) - mu(sp2)).^2 ./ sig(sp2).^2);
end
semilogx(r_arr, mu,'k--', 'linewidth',3); hold on
xlim([4e-1,1e3])
xlabel('arcsec', 'fontsize',15);
ylabel('<I>_{stack}','fontsize',15);
title('CIBER galaxies','fontsize', 15)
subplot(2,4,7)
plot(chi0_arr./numel(r_arr),'k.','markersize',20);hold on
plot(chi1_arr./numel(sp1),'r.','markersize',20);
plot(chi2_arr./numel(sp2),'b.','markersize',20);
xlim([0.8,10.2])
ylabel('\chi^2 / N', 'fontsize', 15)
xlabel('sub-stack', 'fontsize', 15)


subplot(2,4,4)
mu = mean(dat.profgps_mat);
sig = nanstd(dat.profgps_mat);
chi0_arr = zeros([1, dat.Nsubs]);
chi1_arr = zeros([1, dat.Nsubs]);
chi2_arr = zeros([1, dat.Nsubs]);
for i = 1:dat.Nsubs
    semilogx(r_arr, dat.profgps_mat(i,:)); hold on
    chi0_arr(i) = sum((dat.profgps_mat(i,:) - mu).^2 ./ sig.^2);
    chi1_arr(i) = sum((dat.profgps_mat(i,sp1) - mu(sp1)).^2 ./ sig(sp1).^2);
    chi2_arr(i) = sum((dat.profgps_mat(i,sp2) - mu(sp2)).^2 ./ sig(sp2).^2);
end
semilogx(r_arr, mu,'k--', 'linewidth',3); hold on
xlim([4e-1,1e3])
xlabel('arcsec', 'fontsize',15);
ylabel('<I>_{stack}','fontsize',15);
title('Sim galaxies','fontsize', 15)
subplot(2,4,8)
plot(chi0_arr./numel(r_arr),'k.','markersize',20);hold on
plot(chi1_arr./numel(sp1),'r.','markersize',20);
plot(chi2_arr./numel(sp2),'b.','markersize',20);
xlim([0.8,10.2])
ylabel('\chi^2 / N', 'fontsize', 15)
xlabel('sub-stack', 'fontsize', 15)

legend_name = {};
legend_name{end+1} = sprintf('all (%d bins)',numel(r_arr));
legend_name{end+1} = sprintf('r(mask) < r < 100\''\'' (%d bins)',numel(sp1));
legend_name{end+1} = sprintf('r > 100\''\'' (%d bins)',numel(sp2));
legend(legend_name);

savename=strcat(pltsavedir,'chi_sub_profile');
%print(savename,'-dpng');%close
%%

cbexcess_mat = dat.profgcb_mat - dat.profscb_mat;
psexcess_mat = dat.profgps_mat - dat.profsps_mat;
excess_mat = cbexcess_mat - psexcess_mat;

figure
setwinsize(gcf,900,400)

subplot(2,3,1)
mu = mean(cbexcess_mat);
sig = nanstd(cbexcess_mat);
chi0_arr = zeros([1, dat.Nsubs]);
chi1_arr = zeros([1, dat.Nsubs]);
chi2_arr = zeros([1, dat.Nsubs]);
for i = 1:dat.Nsubs
    semilogx(r_arr, cbexcess_mat(i,:)); hold on
    chi0_arr(i) = sum((cbexcess_mat(i,:) - mu).^2 ./ sig.^2);
    chi1_arr(i) = sum((cbexcess_mat(i,sp1) - mu(sp1)).^2 ./ sig(sp1).^2);
    chi2_arr(i) = sum((cbexcess_mat(i,sp2) - mu(sp2)).^2 ./ sig(sp2).^2);
end
semilogx(r_arr, mu,'k--', 'linewidth',3); hold on
xlim([4e-1,1e3])
xlabel('arcsec', 'fontsize',15);
ylabel('<I>_{stack}','fontsize',15);
title('CIBER excess','fontsize', 15);
%legend('show');
subplot(2,3,4)
plot(chi0_arr./numel(r_arr),'k.','markersize',20);hold on
plot(chi1_arr./numel(sp1),'r.','markersize',20);
plot(chi2_arr./numel(sp2),'b.','markersize',20);
xlim([0.8,10.2])
ylabel('\chi^2 / N', 'fontsize', 15)
xlabel('sub-stack', 'fontsize', 15)

subplot(2,3,2)
mu = mean(psexcess_mat);
sig = nanstd(psexcess_mat);
chi0_arr = zeros([1, dat.Nsubs]);
chi1_arr = zeros([1, dat.Nsubs]);
chi2_arr = zeros([1, dat.Nsubs]);
for i = 1:dat.Nsubs
    semilogx(r_arr, psexcess_mat(i,:)); hold on
    chi0_arr(i) = sum((psexcess_mat(i,:) - mu).^2 ./ sig.^2);
    chi1_arr(i) = sum((psexcess_mat(i,sp1) - mu(sp1)).^2 ./ sig(sp1).^2);
    chi2_arr(i) = sum((psexcess_mat(i,sp2) - mu(sp2)).^2 ./ sig(sp2).^2);
end
semilogx(r_arr, mu,'k--', 'linewidth',3); hold on
xlim([4e-1,1e3])
xlabel('arcsec', 'fontsize',15);
ylabel('<I>_{stack}','fontsize',15);
title('Sim excess','fontsize', 15);
%legend('show');
subplot(2,3,5)
plot(chi0_arr./numel(r_arr),'k.','markersize',20);hold on
plot(chi1_arr./numel(sp1),'r.','markersize',20);
plot(chi2_arr./numel(sp2),'b.','markersize',20);
xlim([0.8,10.2])
ylabel('\chi^2 / N', 'fontsize', 15)
xlabel('sub-stack', 'fontsize', 15)

subplot(2,3,3)
mu = mean(excess_mat);
sig = nanstd(excess_mat);
chi0_arr = zeros([1, dat.Nsubs]);
chi1_arr = zeros([1, dat.Nsubs]);
chi2_arr = zeros([1, dat.Nsubs]);
for i = 1:dat.Nsubs
    semilogx(r_arr, excess_mat(i,:)); hold on
    chi0_arr(i) = sum((excess_mat(i,:) - mu).^2 ./ sig.^2);
    chi1_arr(i) = sum((excess_mat(i,sp1) - mu(sp1)).^2 ./ sig(sp1).^2);
    chi2_arr(i) = sum((excess_mat(i,sp2) - mu(sp2)).^2 ./ sig(sp2).^2);
end
semilogx(r_arr, mu,'k--', 'linewidth',3); hold on
xlim([4e-1,1e3])
xlabel('arcsec', 'fontsize',15);
ylabel('<I>_{stack}','fontsize',15);
title('Excess','fontsize', 15);
%legend('show');
subplot(2,3,6)
plot(chi0_arr./numel(r_arr),'k.','markersize',20);hold on
plot(chi1_arr./numel(sp1),'r.','markersize',20);
plot(chi2_arr./numel(sp2),'b.','markersize',20);
xlim([0.8,10.2])
ylabel('\chi^2 / N', 'fontsize', 15)
xlabel('sub-stack', 'fontsize', 15)

savename=strcat(pltsavedir,'chi_sub_excess');
print(savename,'-dpng');%close
