flight = 40030;
inst = 1;

mypaths=get_paths(flight);
figure
setwinsize(gcf,1200,600)

for im = 1:3

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

r_arr = ihlprofdat.r_arr;
m_min = ihlprofdat.m_min_arr(im);
m_max = ihlprofdat.m_max_arr(im);
    
rmask = 0;
for i = 4:8
    rmask = rmask + get_mask_radius(inst,i,(m_min + m_max)/2);
end
rmask = rmask ./ 5;
sp1 = find(r_arr > rmask & r_arr <= 100);
sp2 = find(r_arr > 100);

chi1dat_arr = zeros([3,5]);
chi2dat_arr = zeros([3,5]);
name = {};
for ifield = 4:8
    dt=get_dark_times(flight,inst,ifield);
    name{end+1} = dt.name;
    loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
    load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');    
    diff = ihlprofdat.excess(im).diff;
    diff_err = ihlprofdat.excess(im).diff_err;
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
legend(legend_name);
set(gca,'xTick',4:8); 
set(gca, 'xTickLabels', name); 
xlim([3.8,8.2])
title(strcat(num2str(m_min),' < m < ',num2str(m_max)),'fontsize', 15)
ylabel('\chi / N', 'fontsize', 15)
subplot(2,3,im+3)
semilogy(4:8,chi2dat_arr(1,:)./numel(r_arr),'k.','markersize',20);hold on
plot(4:8,chi2dat_arr(2,:)./numel(sp1),'r.','markersize',20);
plot(4:8,chi2dat_arr(3,:)./numel(sp2),'b.','markersize',20);
set(gca,'xTick',4:8);
set(gca, 'xTickLabels', name); 
xlim([3.8,8.2])
title(strcat(num2str(m_min),' < m < ',num2str(m_max)),'fontsize', 15)
ylabel('\chi^2 / N', 'fontsize', 15)
end
pltsavedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));
savename=strcat(pltsavedir,'excess_chi');
%print(savename,'-dpng');%close
