flight = 40030;
inst = 1;
mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

nbins = 25;

figure
setwinsize(gcf,1200,600)
idx_set = nchoosek(4:8,2);
t_mat = zeros([nbins,size(idx_set,1)]);
for im = 1:3
xtick_name = {};
for i = 1:size(idx_set,1)
    idx1 = idx_set(i,1);
    idx2 = idx_set(i,2);
    dt1=get_dark_times(flight,inst,idx1);
    loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
    load(sprintf('%s/%s_ihlprofdat',loaddir,dt1.name),'ihlprofdat');
    diff1 = ihlprofdat.excess(im).diff;
    diff_err1 = ihlprofdat.excess(im).diff_err;
    dt2=get_dark_times(flight,inst,idx2);
    loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
    load(sprintf('%s/%s_ihlprofdat',loaddir,dt2.name),'ihlprofdat');
    diff2 = ihlprofdat.excess(im).diff;
    diff_err2 = ihlprofdat.excess(im).diff_err;

    t = abs((diff1 - diff2)./sqrt(diff_err1.^2 + diff_err2.^2));
    t_mat(:,i) = t;
    
    xtick_name{end+1} = strcat(dt1.name,'--',dt2.name);
end

subplot(2,3,im)
legend_name = {};
for ibin = 2:8:nbins
    plot(1:size(idx_set,1),t_mat(ibin,:),'linewidth',2);hold on
    legend_name{end+1} = sprintf('r = %.2f\''\''',ihlprofdat.r_arr(ibin));
end
legend(legend_name);
xlim([0.5,10.5])
set(gca,'xTick',1:size(idx_set,1)); 
set(gca, 'xTickLabels', xtick_name); 
xtickangle(45)
ylabel('t-value', 'fontsize', 15)

subplot(2,3,im+3)
legend_name = {};
for ibin = 2:nbins
    plot(1:size(idx_set,1),t_mat(ibin,:),'k','linewidth',2);hold on
end
xlim([0.5,10.5])
set(gca,'xTick',1:size(idx_set,1)); 
set(gca, 'xTickLabels', xtick_name); 
xtickangle(45)
ylabel('t-value', 'fontsize', 15)
end

savename=strcat(pltsavedir,'t_test');
print(savename,'-dpng');%close
