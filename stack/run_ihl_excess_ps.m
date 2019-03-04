flight=40030;
inst=1;
ifield=4;

dt=get_dark_times(flight,inst,ifield);
loaddir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/');
savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

load(sprintf('%s/ciber_ps/%s_profstack',loaddir,dt.name));

m_min_arr = [0,8:22];
m_max_arr = [8:23];

% fit a power law between 50-500 arcsec
fitmin = 8;
fitmax = 50;

for im=10:13%1:numel(m_min_arr)
    m_min = profstack(im).m_min;
    m_max = profstack(im).m_max;
    r_arr = profstack(im).r_arr;
    Ns = profstack(im).Ns;
    Ng = profstack(im).Ng;
    diffcb = profstack(im).profgcb - profstack(im).profscb;
    diffcb_err = sqrt(profstack(im).profgcb_err.^2 + profstack(im).profscb_err.^2);
    diffps = profstack(im).profgps - profstack(im).profsps;
    diffps_err = sqrt(profstack(im).profgps_err.^2 + profstack(im).profsps_err.^2);
    
    diff = diffcb - diffps;
    diff_err = sqrt(diffcb_err.^2 + diffps_err.^2);
        
    [fit,fit_model] = fit_pl_prof_chi2(r_arr,diff,diff_err,fitmin,fitmax);
    
    profstack(im).diff = diff;
    profstack(im).diff_err = diff_err;
    profstack(im).pl_fit_params = fit.params;
    profstack(im).pl_fit = fit.params(1).*r_arr.^fit.params(2);
    

    figure
    setwinsize(gcf,1000,400)
    
    subplot(1,2,1)
    errorbar(r_arr.*0.98,diffcb,diffcb_err,'r.','markersize',10,...
        'DisplayName','CB gals-stars');hold on
    errorbar(r_arr.*1.02,diffps,diffps_err,'b.','markersize',10,...
        'DisplayName','PS gals-stars');
    errorbar(r_arr,diff,diff_err,'ko','DisplayName',...
        'CIBER excess: (CB gals-stars) - (PS gals-stars)');
    ylim1 = gca; ylim1 = ylim1.YLim;
    plot(r_arr,fit_model,'k--','DisplayName',...
        strcat('power low fit (',num2str(fitmin),'<r<',num2str(fitmax),' arcsec)'));
    
    subplot(1,2,2)
    loglog(r_arr.*0.98,diffcb,'r.','markersize',10);hold on
    loglog(r_arr.*1.02,diffps,'b.','markersize',10);
    loglog(r_arr,diff,'k.','markersize',10);
    errorbar(r_arr.*0.98,diffcb,diffcb_err,'r.');
    errorbar(r_arr.*1.02,diffps,diffps_err,'b.');
    errorbar(r_arr,diff,diff_err,'k.');
    
    loglog(r_arr.*0.98,-diffcb,'ro','markersize',5);hold on
    loglog(r_arr.*1.02,-diffps,'bo','markersize',5);
    loglog(r_arr,-diff,'ko','markersize',5);
    errorbar(r_arr.*0.98,-diffcb,diffcb_err,'ro','markersize',5);
    errorbar(r_arr.*1.02,-diffps,diffps_err,'bo','markersize',5);
    errorbar(r_arr,-diff,diff_err,'ko','markersize',5);
    ylim2 = gca; ylim2 = ylim2.YLim;
    loglog(r_arr,fit_model,'k--');

    subplot(1,2,1)
    xlabel('arcsec')
    ylabel('normalized <I_{stack}>')
    xlim([4e-1,1e3])
    ylim(ylim1)
    set(gca, 'XScale', 'log')
    h=legend('show','Location','southeast');
    set(h,'fontsize',10)
    legend boxoff

    subplot(1,2,2)
    xlim([4e-1,1e3])
    ylim(ylim2)
    xlabel('arcsec')
    ylabel('normalized <I_{stack}>')
    
    suptitle(strcat(num2str(m_min),'<mAB(y band)<',num2str(m_max)));
    
    savename=strcat(savedir,dt.name,'_rprof',num2str(m_min),'_',num2str(m_max),'diff');
    print(savename,'-dpng');%close

end
save(sprintf('%s/ciber_ps/%s_profstack',loaddir,dt.name),'profstack');
%% compare diff for all fields

flight=40030;
inst=1;

loaddir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/');
savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

m_min_arr = [0,8:22];
m_max_arr = [8:23];


h=zeros(1,5);
for im=10:13
m_min = m_min_arr(im);
m_max = m_max_arr(im);
figure
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    load(sprintf('%s/ciber_ps/%s_profstack',loaddir,dt.name));
    r_arr = profstack(im).r_arr;
    diff = profstack(im).diff;
    diff_err = profstack(im).diff_err;
    h(ifield-3) = semilogx(r_arr,diff,'.','color',get_color(ifield-3),...
        'markersize',10,'DisplayName',dt.name);hold on
    errorbar(r_arr,diff,diff_err,'.','color',get_color(ifield-3));
end
ylim([-0.2,0.2])
xlim([4e-1,1e3])
leg=legend(h,'Location','southeast');
set(leg,'fontsize',10)
legend boxoff
title(strcat(num2str(m_min),'<mAB(y band)<',num2str(m_max)));

savename=strcat(savedir,'allfield_rprof',num2str(m_min),'_',num2str(m_max),'diff');
print(savename,'-dpng');%close

end