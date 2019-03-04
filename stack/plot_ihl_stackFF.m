function plot_ihl_stackFF(flight,inst,ifield)

mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);


pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdatFF',loaddir,dt.name),'ihlprofdat');

%%% write the basic info %%%
r_arr = ihlprofdat.r_arr;
psf_arr = ihlprofdat.psf_arr;
m_min_arr = ihlprofdat.m_min_arr;
m_max_arr = ihlprofdat.m_max_arr;
counts_arr = ihlprofdat.counts_arr;
countg_arr = ihlprofdat.countg_arr;

%%% plot the stacking profile %%%
figure
setwinsize(gcf,1000,600)

for im=10:12
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    counts = counts_arr(im);
    countg = countg_arr(im);
    
    profscb = ihlprofdat.norm(im).profscb;
    profscb_err = ihlprofdat.norm(im).profscb_err;
    profsps = ihlprofdat.norm(im).profsps;
    profsps_err = ihlprofdat.norm(im).profsps_err;
    profgcb = ihlprofdat.norm(im).profgcb;
    profgcb_err = ihlprofdat.norm(im).profgcb_err;
    profgps = ihlprofdat.norm(im).profgps;
    profgps_err = ihlprofdat.norm(im).profgps_err;

    subplot(2,3,im-9)
    loglog(r_arr.*0.98,profscb,'r.','markersize',10,'DisplayName','stars w/ err');
    hold on
    loglog(r_arr.*1.02,profgcb,'b.','markersize',10,'DisplayName','gals w/ err');
    loglog(r_arr.*0.98,profsps,'m.','markersize',10,'DisplayName','stars');
    loglog(r_arr.*1.02,profgps,'c.','markersize',10,'DisplayName','gals');
    h=legend('show','Location','northeast');
    set(h,'fontsize',10)
    legend boxoff

    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    errorbar(r_arr.*0.98,profsps,profsps_err,'m.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'c.','markersize',10);    
    loglog(r_arr.*0.98,-profscb,'ro','markersize',5);hold on
    loglog(r_arr.*1.02,-profgcb,'bo','markersize',5);
    loglog(r_arr.*0.98,-profsps,'mo','markersize',5);hold on
    loglog(r_arr.*1.02,-profgps,'co','markersize',5);
    errorbar(r_arr.*0.98,-profscb,profscb_err,'ro','markersize',5);
    errorbar(r_arr.*1.02,-profgcb,profgcb_err,'bo','markersize',5);
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
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
        ' (',num2str(counts),' stars, ',num2str(countg), ' gals)'),'fontsize',15)
    
    subplot(2,3,im-6)
    diffcb = ihlprofdat.excess(im).diffcb;
    diffcb_err = ihlprofdat.excess(im).diffcb_err;
    diffps = ihlprofdat.excess(im).diffps;
    diffps_err = ihlprofdat.excess(im).diffsp_err;
    diff = ihlprofdat.excess(im).diff;
    diff_err = ihlprofdat.excess(im).diff_err;
    
    loglog(r_arr.*0.98,diffcb,'r.','markersize',10,...
        'DisplayName','gals - stars (w/ FF err)');hold on
    loglog(r_arr.*1.02,diffps,'b.','markersize',10,...
        'DisplayName','gals - stars');
    h=legend('show','Location','northeast');
    set(h,'fontsize',10)
    legend boxoff    
    errorbar(r_arr.*0.98,diffcb,diffcb_err,'r.');
    errorbar(r_arr.*1.02,diffps,diffps_err,'b.');
    

    loglog(r_arr.*0.98,-diffcb,'ro','markersize',5);hold on
    loglog(r_arr.*1.02,-diffps,'bo','markersize',5);
    errorbar(r_arr.*0.98,-diffcb,diffcb_err,'ro','markersize',5);
    errorbar(r_arr.*1.02,-diffps,diffps_err,'bo','markersize',5);
    
    xlim([4e-1,1e3])
    ylim([1e-6,1e0])
    xlabel('arcsec')
    ylabel('<I_{stack}>')
end

    savename=strcat(pltsavedir,dt.name,'_rprof_excessFF');
    print(savename,'-dpng');close

%{
for im=10:12
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    counts = counts_arr(im);
    countg = countg_arr(im);

    figure
    setwinsize(gcf,1000,800)    
    profscb = ihlprofdat.data(im).profscb;
    profscb_err = ihlprofdat.data(im).profscb_err;
    profsps = ihlprofdat.data(im).profsps;
    profsps_err = ihlprofdat.data(im).profsps_err;
    profgcb = ihlprofdat.data(im).profgcb;
    profgcb_err = ihlprofdat.data(im).profgcb_err;
    profgps = ihlprofdat.data(im).profgps;
    profgps_err = ihlprofdat.data(im).profgps_err;
    
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

    profscb = ihlprofdat.norm(im).profscb;
    profscb_err = ihlprofdat.norm(im).profscb_err;
    profsps = ihlprofdat.norm(im).profsps;
    profsps_err = ihlprofdat.norm(im).profsps_err;
    profgcb = ihlprofdat.norm(im).profgcb;
    profgcb_err = ihlprofdat.norm(im).profgcb_err;
    profgps = ihlprofdat.norm(im).profgps;
    profgps_err = ihlprofdat.norm(im).profgps_err;
    
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
    
%     savename=strcat(pltsavedir,...
%         dt.name,'_rprof',num2str(m_min),'_',num2str(m_max),'bksub');
%     print(savename,'-dpng');close
end
%}

%%% plot the excess profile %%%
%{
figure
setwinsize(gcf,1200,300)

for im=10:12
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);

    diffcb = ihlprofdat.excess(im).diffcb;
    diffcb_err = ihlprofdat.excess(im).diffcb_err;
    diffps = ihlprofdat.excess(im).diffps;
    diffps_err = ihlprofdat.excess(im).diffsp_err;
    diff = ihlprofdat.excess(im).diff;
    diff_err = ihlprofdat.excess(im).diff_err;
    
    subplot(1,3,im-9)
    loglog(r_arr.*0.98,diffcb,'r.','markersize',10,...
        'DisplayName','CB diff = CB gal - CB stars');hold on
    loglog(r_arr.*1.02,diffps,'b.','markersize',10,...
        'DisplayName','PS diff = PS gal - PS stars');
    %loglog(r_arr,diff,'k.','markersize',10,...
    %    'DisplayName','CB diff - PS diff');
    h=legend('show','Location','southeast');
    set(h,'fontsize',10)
    legend boxoff    
    errorbar(r_arr.*0.98,diffcb,diffcb_err,'r.');
    errorbar(r_arr.*1.02,diffps,diffps_err,'b.');
    %errorbar(r_arr,diff,diff_err,'k.');
    

    loglog(r_arr.*0.98,-diffcb,'ro','markersize',5);hold on
    loglog(r_arr.*1.02,-diffps,'bo','markersize',5);
    %loglog(r_arr,-diff,'ko','markersize',5);
    errorbar(r_arr.*0.98,-diffcb,diffcb_err,'ro','markersize',5);
    errorbar(r_arr.*1.02,-diffps,diffps_err,'bo','markersize',5);
    %errorbar(r_arr,-diff,diff_err,'ko','markersize',5);
    
    xlim([4e-1,1e3])
    ylim([1e-6,1e0])
    xlabel('arcsec')
    ylabel('normalized <I_{stack}>')
    
    title(strcat(num2str(m_min),'<mAB(y band)<',num2str(m_max)));
end
% savename=strcat(pltsavedir,dt.name,'_rprof_excess');
% print(savename,'-dpng');close
%}

return