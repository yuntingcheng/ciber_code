function plot_ihl_stack_hist_sub(flight,inst)

mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

%%% plot the excess profile in each field %%%
%%
for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);
figure
setwinsize(gcf,1300,1200)

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');
r_arr = ihlprofdat.r_arr;
ihlall = zeros([3,numel(r_arr)]);
for im=1:3
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1; 
    ihlall(im,:) = ihl;
    subplot(4,3,im)
    semilogx(r_arr,ihl./ihlall(im,:),'k.','markersize',10);hold on
    errorbar(r_arr,ihl./ihlall(im,:),ihl_err./ihlall(im,:),'k.','markersize',10);
    ylim([1-nanmedian(ihl_err)*5,1+nanmedian(ihl_err)*5]);
    if im==2
        title('all (k), middle (r), edge (b)');

    end
    subplot(4,3,im+3)
    semilogx(r_arr,ihl./ihlall(im,:),'k.','markersize',10);hold on
    errorbar(r_arr,ihl./ihlall(im,:),ihl_err./ihlall(im,:),'k.','markersize',10); 
    ylim([1-nanmedian(ihl_err)*5,1+nanmedian(ihl_err)*5]);

    if im==2
        title('all (k), top (r), bottom (b)');
    end
    subplot(4,3,im+6)
    semilogx(r_arr,ihl./ihlall(im,:),'k.','markersize',10);hold on
    errorbar(r_arr,ihl./ihlall(im,:),ihl_err./ihlall(im,:),'k.','markersize',10);  
    ylim([1-nanmedian(ihl_err)*5,1+nanmedian(ihl_err)*5]);
    if im==2
        title('all (k), left (r), right (b)');
    end
    subplot(4,3,im+9)
    semilogx(r_arr,ihl./ihlall(im,:),'k.','markersize',10);hold on
    errorbar(r_arr,ihl./ihlall(im,:),ihl_err./ihlall(im,:),'k.','markersize',10); 
    ylim([1-nanmedian(ihl_err)*5,1+nanmedian(ihl_err)*5]);
    if im==2
        title('all (k), Q1 (r), Q2 (b), Q3 (m), Q4 (c)');
    end
end
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat_hist_%s',loaddir,dt.name,'M'),'ihlprofdat');
for im=1:3
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;    
    subplot(4,3,im)
    h(2)=semilogx(r_arr.*1.05,ihl./ihlall(im,:),'r.','markersize',10);hold on
    errorbar(r_arr.*1.05,ihl./ihlall(im,:),...
        ihl_err./ihlall(im,:),'r.','markersize',10); 
end
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat_hist_%s',loaddir,dt.name,'E'),'ihlprofdat');
for im=1:3
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;
    subplot(4,3,im)
    semilogx(r_arr.*0.95,ihl./ihlall(im,:),'b.','markersize',10);hold on
    errorbar(r_arr.*0.95,ihl./ihlall(im,:),...
        ihl_err./ihlall(im,:),'b.','markersize',10);    
    xlim([7e-1,1.3e3])
    ylabel('I / I (all)', 'fontsize',15)
end

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat_hist_%s',loaddir,dt.name,'T'),'ihlprofdat');
for im=1:3
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;    
    subplot(4,3,im+3)
    semilogx(r_arr.*1.05,ihl./ihlall(im,:),'r.','markersize',10);hold on
    errorbar(r_arr.*1.05,ihl./ihlall(im,:),...
        ihl_err./ihlall(im,:),'r.','markersize',10);    
end
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat_hist_%s',loaddir,dt.name,'B'),'ihlprofdat');
for im=1:3
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;
    subplot(4,3,im+3)
    semilogx(r_arr.*0.95,ihl./ihlall(im,:),'b.','markersize',10);hold on
    errorbar(r_arr.*0.95,ihl./ihlall(im,:),...
        ihl_err./ihlall(im,:),'b.','markersize',10);    
    xlim([7e-1,1.3e3])
    ylabel('I / I (all)', 'fontsize',15)
end

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat_hist_%s',loaddir,dt.name,'L'),'ihlprofdat');
for im=1:3
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;    
    subplot(4,3,im+6)
    semilogx(r_arr.*1.05,ihl./ihlall(im,:),'r.','markersize',10);hold on
    errorbar(r_arr.*1.05,ihl./ihlall(im,:),...
        ihl_err./ihlall(im,:),'r.','markersize',10);    
end
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat_hist_%s',loaddir,dt.name,'R'),'ihlprofdat');
for im=1:3
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;
    subplot(4,3,im+6)
    semilogx(r_arr.*0.95,ihl./ihlall(im,:),'b.','markersize',10);hold on
    errorbar(r_arr.*0.95,ihl./ihlall(im,:),...
        ihl_err./ihlall(im,:),'b.','markersize',10);    
    xlim([7e-1,1.3e3])
    ylabel('I / I (all)', 'fontsize',15)
end

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat_hist_%s',loaddir,dt.name,'Q1'),'ihlprofdat');
for im=1:3
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;    
    subplot(4,3,im+9)
    semilogx(r_arr.*1.05,ihl./ihlall(im,:),'r.','markersize',10);hold on
    errorbar(r_arr.*1.05,ihl./ihlall(im,:),...
        ihl_err./ihlall(im,:),'r.','markersize',10);    
end
load(sprintf('%s/%s_ihlprofdat_hist_%s',loaddir,dt.name,'Q2'),'ihlprofdat');
for im=1:3
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;    
    subplot(4,3,im+9)
    semilogx(r_arr.*1.02,ihl./ihlall(im,:),'b.','markersize',10);hold on
    errorbar(r_arr.*1.02,ihl./ihlall(im,:),...
        ihl_err./ihlall(im,:),'b.','markersize',10);    
end
load(sprintf('%s/%s_ihlprofdat_hist_%s',loaddir,dt.name,'Q3'),'ihlprofdat');
for im=1:3
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;    
    subplot(4,3,im+9)
    h(4)=semilogx(r_arr.*0.98,ihl./ihlall(im,:),'m.','markersize',10);hold on
    errorbar(r_arr.*0.98,ihl./ihlall(im,:),...
        ihl_err./ihlall(im,:),'m.','markersize',10);    
end
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat_hist_%s',loaddir,dt.name,'Q4'),'ihlprofdat');
for im=1:3
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;
    subplot(4,3,im+9)
    h(5)=semilogx(r_arr.*0.95,ihl./ihlall(im,:),'c.','markersize',10);hold on
    errorbar(r_arr.*0.95,ihl./ihlall(im,:),...
        ihl_err./ihlall(im,:),'c.','markersize',10);    
    xlim([7e-1,1.3e3])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I / I (all)', 'fontsize',15)
end

suptitle(dt.name);
savename=strcat(pltsavedir,dt.name,'_excessprof_hist_subreg');
print(savename,'-dpng');%close

end
%%
return