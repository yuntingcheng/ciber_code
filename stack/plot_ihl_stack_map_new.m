flight=40030;
inst=1;
masklim = false;
savefig = false;
mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.alldat,'plots/TM',num2str(inst),'/'));
%% jackknife prof, Src prof, BG prof
for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if masklim
    load(sprintf('%s/stackdat_%s_masklim',...
            loaddir,dt.name),'stackdatall');        
    load(sprintf('%s/excessdat_%s_masklim',...
            loaddir,dt.name),'excessdatall');        
else
    load(sprintf('%s/stackdat_%s',...
            loaddir,dt.name),'stackdatall');
    load(sprintf('%s/excessdat_%s',...
            loaddir,dt.name),'excessdatall');
end

figure
setwinsize(gcf,1200,500)
for im=1:4
    stackdat = stackdatall(im).stackdat;
    excessdat = excessdatall(im).excessdat;
    m_min = stackdat.m_min;
    m_max = stackdat.m_max;
    countg = stackdat.all.countg;
    r_arr = stackdat.r_arr;
    
    %%% jackknife, Src, BG inner %%%
    subplot(2,4,im)
    semilogx(r_arr,stackdat.jack(1).profcbg,'color',[0.8,1,0.8]);hold on
    
    prof = stackdat.all.profcbg;
    err = sqrt(diag(excessdat.datcov.covcb));
    errorbar(r_arr, prof, err, 'b.','markersize',10); 
    profbg = stackdat.bg.profcbg;
    bg_err = sqrt(diag(excessdat.bgcov.covcb));
    errorbar(r_arr.*1.03, profbg, bg_err, 'k.-','markersize',5);
    h=legend({'jackknife profile','Src profile','BG profile'},...
        'Location','northeast');
    set(h,'fontsize',7)
    for i=2:numel(stackdat.jack)
        semilogx(r_arr,stackdat.jack(i).profcbg,'color',[0.4,1,1]);
    end
    errorbar(r_arr, prof, err, 'b.','markersize',10);
    errorbar(r_arr.*1.03, profbg, bg_err, 'k.','markersize',10);
    xlim([4e-1,50])
    if im==1
        ylabel('I [nW/m^2/sr]', 'fontsize',15)
    end
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
        ' (',num2str(countg), ' galaxies)'),'fontsize',15)
    
    %%% jackknife, Src, BG outer %%%
    subplot(2,4,im+4)
    semilogx(r_arr,stackdat.jack(1).profcbg,'color',[0.8,1,0.8]);hold on
    
    prof = stackdat.all.profcbg;
    err = sqrt(diag(excessdat.datcov.covcb));
    profbg = stackdat.bg.profcbg;
    bg_err = sqrt(diag(excessdat.bgcov.covcb));
    for i=2:numel(stackdat.jack)
        semilogx(r_arr,stackdat.jack(i).profcbg,'color',[0.4,1,1]);
    end
    errorbar(r_arr, prof, err, 'b.','markersize',10);
    errorbar(r_arr.*1.03, profbg, bg_err, 'k.','markersize',10);
    xlim([20,1.1e3])
    if im==1
        ylabel('I [nW/m^2/sr]', 'fontsize',15)
    end
    xlabel('r [arcsec]', 'fontsize',15)  
end

savename = sprintf('%s%s_profs',pltsavedir,dt.name);
if masklim
    suptitle(strcat(dt.name,' (mask to mag bin max)'));
    savename = strcat(savename,'_masklim'); 
else
    suptitle(strcat(dt.name,' (mask all PanSTARRS sources)'));
end
if savefig
    print(savename,'-dpng');close
end
end
%% BG-sub prof, scaled PSF

for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if masklim
    load(sprintf('%s/stackdat_%s_masklim',...
            loaddir,dt.name),'stackdatall');        
    load(sprintf('%s/excessdat_%s_masklim',...
            loaddir,dt.name),'excessdatall');        
else
    load(sprintf('%s/stackdat_%s',...
            loaddir,dt.name),'stackdatall');
    load(sprintf('%s/excessdat_%s',...
            loaddir,dt.name),'excessdatall');
end
load(sprintf('%s/psfdat_%s',loaddir,dt.name),'psfdatall');
psfdat = psfdatall.comb(im);

figure
setwinsize(gcf,1200,500)
for im=1:4
    stackdat = stackdatall(im).stackdat;
    excessdat = excessdatall(im).excessdat;
   
    m_min = stackdat.m_min;
    m_max = stackdat.m_max;
    countg = stackdat.all.countg;
    r_arr = stackdat.r_arr;
    
    %%% jackknife, Src, BG inner %%%
    subplot(2,4,im)
    prof = stackdat.all.profcbg-stackdat.all.profcbgbg;
    err = sqrt(diag(excessdat.datcov.covcb+excessdat.bgcov.covcb));
    semilogx(r_arr, prof,'b.-','markersize',10);hold on
    
    profpsf = psfdat.all.profcb.*excessdat.normcb;
    errpsf = sqrt(diag(excessdat.psfcov.covcb));
    errorbar(r_arr.*1.02, profpsf, errpsf, 'r','markersize',10);
    
    h=legend({'BG-sub Src profile','scaled PSF'},...
        'Location','northeast');
    set(h,'fontsize',7)

    errorbar(r_arr, prof, err, 'b.','markersize',10);
    xlim([4e-1,50])
    if im==1
        ylabel('I [nW/m^2/sr]', 'fontsize',15)
    end
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
        ' (',num2str(countg), ' galaxies)'),'fontsize',15)
    
    subplot(2,4,im+4)
    prof = stackdat.all.profcbg-stackdat.all.profcbgbg;
    err = sqrt(diag(excessdat.datcov.covcb+excessdat.bgcov.covcb));
    semilogx(r_arr, prof,'b.-','markersize',10);hold on
    profpsf = psfdat.all.profcb.*excessdat.normcb;
    errpsf = sqrt(diag(excessdat.psfcov.covcb));
    errorbar(r_arr.*1.02, profpsf, errpsf, 'r','markersize',10);
    errorbar(r_arr, prof, err, 'b.','markersize',10);
    xlim([10,1.1e3])
    if im==1
        ylabel('I [nW/m^2/sr]', 'fontsize',15)
    end
    xlabel('r [arcsec]', 'fontsize',15)  
    
end

savename = sprintf('%s%s_bgsubprofs',pltsavedir,dt.name);
if masklim
    suptitle(strcat(dt.name,' (mask to mag bin max)'));
    savename = strcat(savename,'_masklim'); 
else
    suptitle(strcat(dt.name,' (mask all PanSTARRS sources)'));
end
if savefig
    print(savename,'-dpng');close
end

end
%% Excess
for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if masklim
    load(sprintf('%s/stackdat_%s_masklim',...
            loaddir,dt.name),'stackdatall');        
    load(sprintf('%s/excessdat_%s_masklim',...
            loaddir,dt.name),'excessdatall');        
else
    load(sprintf('%s/stackdat_%s',...
            loaddir,dt.name),'stackdatall');
    load(sprintf('%s/excessdat_%s',...
            loaddir,dt.name),'excessdatall');
end

figure
setwinsize(gcf,1200,600)
for im=1:4
    stackdat = stackdatall(im).stackdat;
    excessdat = excessdatall(im).excessdat;
    weight = stackdat.all.profhitg;
    hscclusdat = get_hsc_clus_prof(flight, inst, masklim, weight);

    m_min = stackdat.m_min;
    m_max = stackdat.m_max;
    countg = stackdat.all.countg;
    r_arr = stackdat.r_arr;
    
    %%% jackknife, Src, BG inner %%%
    subplot(2,4,im)
    if im==2
        loglog(1e-8,1e-8,'k.');hold on
        plot(1e-8,1e-8,'b.');
        plot(1e-8,1e-8,'r.');
        h=legend({'CIBER Excess','PanSTARRS Excess','HSC clustering Excess'},...
                'Location','northeast');
        set(h,'fontsize',6)
    end
    prof = excessdat.excess.profcbg;
    err = sqrt(diag(excessdat.excov.covcb))';
    err(err > abs(prof)) = abs(prof(err > abs(prof))) - 1e-10;
    loglog(r_arr, prof,'k.','markersize',10);hold on
    errorbar(r_arr, prof,err,'k.','markersize',10);
    errorbar(r_arr, -prof,err,'ko','markersize',5);
    prof = excessdat.excess.profpsg;
    err = sqrt(diag(excessdat.excov.covps))';
    err(err > abs(prof)) = abs(prof(err > abs(prof))) - 1e-10;
    errorbar(r_arr.*1.03, prof,err,'b.','markersize',10);
    errorbar(r_arr.*1.03, -prof,err,'bo','markersize',5);
    prof = hscclusdat(im).prof;
    err = hscclusdat(im).err;
    err(err > abs(prof)) = abs(prof(err > abs(prof))) - 1e-10;
    errorbar(r_arr.*0.97, prof,err,'r.','markersize',10);
    errorbar(r_arr.*0.97, -prof,err,'ro','markersize',5);
    
    xlim([4e-1,1.1e3])
    ylim([1e-4,1e3])
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
        ' (',num2str(countg), ' galaxies)'),'fontsize',15)
    if im==1
        ylabel('I [nW/m^2/sr]', 'fontsize',15);
    end

    subplot(2,4,im+4)
    prof = excessdat.excess.profcbg;
    err = sqrt(diag(excessdat.excov.covcb))';
    errorbar(r_arr, prof,err,'k.','markersize',10);hold on
    sp = find(r_arr > 100);
    v = excessdat.excess.profcbg100;
    e = excessdat.excov.covcb100;
    plot(r_arr(sp),v*ones(size(sp)),'color',[0.2,0.4,0.2],'linewidth',2);
    if im==1
        h=legend({'CIBER Excess','CIBER Excess > 100'},...
            'Location','northeast');
        set(h,'fontsize',8)
    end
    fill([r_arr(sp),flip(r_arr(sp))],...
    [(v+e) * ones(size(sp)),(v-e) * ones(size(sp))],...
        [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
    set(gca, 'XScale', 'log');
    hline(0,'k:');
    xlim([4e-1,1.1e3])
    ylim([-2,8])
    xlabel('r [arcsec]', 'fontsize',15);
    if im==1
        ylabel('I [nW/m^2/sr]', 'fontsize',15);
    end
       
end

savename = sprintf('%s%s_excess',pltsavedir,dt.name);
if masklim
    suptitle(strcat(dt.name,' (mask to mag bin max)'));
    savename = strcat(savename,'_masklim'); 
else
    suptitle(strcat(dt.name,' (mask all PanSTARRS sources)'));
end
if savefig
    print(savename,'-dpng');close
end

end
%% Cov of individual terms
for ifield = 4%4:8
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if masklim
    load(sprintf('%s/stackdat_%s_masklim',...
            loaddir,dt.name),'stackdatall');        
    load(sprintf('%s/excessdat_%s_masklim',...
            loaddir,dt.name),'excessdatall');        
else
    load(sprintf('%s/stackdat_%s',...
            loaddir,dt.name),'stackdatall');
    load(sprintf('%s/excessdat_%s',...
            loaddir,dt.name),'excessdatall');
end

for itype=1:3
figure
setwinsize(gcf,1400,500)
for im=1:4
    stackdat = stackdatall(im).stackdat;
    excessdat = excessdatall(im).excessdat;
    switch itype
    case 1
        exdat = excessdat.datcov;
        name = 'data';
    case 2
        exdat = excessdat.bgcov;
        name = 'background';
    case 3
        exdat = excessdat.psfcov;
        name = 'PSF';
    end
    
    
    m_min = stackdat.m_min;
    m_max = stackdat.m_max;
    r_arr = stackdat.r_arr;

    subplot(2,4,im)
    imageclip(exdat.covcb);
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);
    xticks([6:6:25])
    xticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
        num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
    xtickangle(45)
    yticks([6:6:25])
    yticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
        num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
    ytickangle(45)
    if im==1
        ylabel(strcat(name, ' Cov'), 'fontsize',15);
    end
    
    subplot(2,4,im+4)
    imageclip(normalize_cov(exdat.covcb));
    xticks([6:6:25])
    xticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
        num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
    xtickangle(45)
    yticks([6:6:25])
    yticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
        num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
    ytickangle(45)
    if im==1
        ylabel(strcat(name, ' Corr'), 'fontsize',15);
    end

end
savename = sprintf('%s%s_cov%s',pltsavedir,dt.name,name);
if masklim
    suptitle(strcat(dt.name,' (mask to mag bin max)'));
    savename = strcat(savename,'_masklim'); 
else
    suptitle(strcat(dt.name,' (mask all PanSTARRS sources)'));
end
if savefig
    print(savename,'-dpng');close
end

end

end

%% Cov old
for ifield = 4%4:8
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if masklim
    load(sprintf('%s/stackdat_%s_masklim',...
            loaddir,dt.name),'stackdatall');        
    load(sprintf('%s/excessdat_%s_masklim',...
            loaddir,dt.name),'excessdatall');        
else
    load(sprintf('%s/stackdat_%s',...
            loaddir,dt.name),'stackdatall');
    load(sprintf('%s/excessdat_%s',...
            loaddir,dt.name),'excessdatall');
end

figure
setwinsize(gcf,1400,750)
for im=1:4
    stackdat = stackdatall(im).stackdat;
    excessdat = excessdatall(im).excessdat;
    m_min = stackdat.m_min;
    m_max = stackdat.m_max;
    r_arr = stackdat.r_arr;

    subplot(3,4,im)
    imageclip(excessdat.datcov.covcb);
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);
    xticks([6:6:25])
    xticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
        num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
    xtickangle(45)
    yticks([6:6:25])
    yticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
        num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
    ytickangle(45)
    if im==1
        ylabel('data Cov', 'fontsize',15);
    end
    
    subplot(3,4,im+4)
    imageclip(excessdat.bgcov.covcb);
    xticks([6:6:25])
    xticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
        num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
    xtickangle(45)
    yticks([6:6:25])
    yticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
        num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
    ytickangle(45)
    if im==1
        ylabel('BG Cov', 'fontsize',15);
    end

    subplot(3,4,im+8)
    dat = excessdat.excess.profcbg;
    edat = sqrt(diag(excessdat.datcov.covcb));
    ebg = sqrt(diag(excessdat.bgcov.covcb));
    epsf = sqrt(diag(excessdat.psfcov.covcb));
    etot = sqrt(diag(excessdat.excov.covcb));
    etotj = sqrt(diag(excessdat.exjcov.covcb));
    etotji = sqrt(diag(excessdat.exjicov.covcb));
    loglog(r_arr, edat,'b--','linewidth',1);hold on
    loglog(r_arr, ebg,'c--','linewidth',1);
    loglog(r_arr, epsf,'r--','linewidth',1);
    loglog(r_arr, etot,'k-','linewidth',2);
    loglog(r_arr, etotj,'k--','linewidth',2);
    loglog(r_arr, etotji,'k:','linewidth',2);
    loglog(r_arr, dat, 'k.', 'markersize', 10);hold on
    if im==4
        h=legend({'data err','BG err','PSF err', 'Excess err', 'Excess Data'},...
                'Location','northeast');
        set(h,'fontsize',6)
    end
    
    loglog(r_arr, -dat, 'ko', 'markersize', 5);hold on
    xlim([4e-1,1.1e3])
    ylim([1e-1,1.1e2])
    xlabel('r [arcsec]', 'fontsize',15);
    if im==1
        ylabel('I [nW/m^2/sr]', 'fontsize',15);
    end

end
savename = sprintf('%s%s_cov',pltsavedir,dt.name);
if masklim
    suptitle(strcat(dt.name,' (mask to mag bin max)'));
    savename = strcat(savename,'_masklim'); 
else
    suptitle(strcat(dt.name,' (mask all PanSTARRS sources)'));
end
if savefig
    print(savename,'-dpng');close
end

end
%% excess cov
for ifield = 4%4:8
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if masklim
    load(sprintf('%s/stackdat_%s_masklim',...
            loaddir,dt.name),'stackdatall');        
else
    load(sprintf('%s/stackdat_%s',...
            loaddir,dt.name),'stackdatall');
end

figure
setwinsize(gcf,1400,750)
for im=1:4
    stackdat = stackdatall(im).stackdat;
    m_min = stackdat.m_min;
    m_max = stackdat.m_max;
    rsub_arr = stackdat.rsub_arr;
    
    subplot(3,4,im)
    imageclip(stackdat.cov.cbsub);
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)
    if im==1
        ylabel('Excess Cov', 'fontsize',15);
    end
    
    subplot(3,4,im+4)
    imageclip(stackdat.cov.cbsub_rho);
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)
    if im==1
        ylabel('Excess Cov (norm.)', 'fontsize',15);
    end
    
    subplot(3,4,im+8)
    imageclip(stackdat.cov.cbsub_inv);
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)
    if im==1
        ylabel('Excess Cov Inv.', 'fontsize',15);
    end
    
end
savename = sprintf('%s%s_excov',pltsavedir,dt.name);
if masklim
    suptitle(strcat(dt.name,' (mask to mag bin max)'));
    savename = strcat(savename,'_masklim'); 
else
    suptitle(strcat(dt.name,' (mask all PanSTARRS sources)'));
end
if savefig
    print(savename,'-dpng');close
end

end

%% example chi2 test
flight=40030;
inst=1;
mypaths=get_paths(flight);
ifield = 4;
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackdat_%s',loaddir,dt.name),'stackdatall');

figure
setwinsize(gcf,1400,300)
im = 2;

stackdat = stackdatall(im).stackdat;
m_min = stackdat.m_min;
m_max = stackdat.m_max;
rsub_arr = stackdat.rsub_arr;

cov = stackdat.cov.cbsub;
covi = stackdat.cov.cbsub_inv;
d = stackdat.excess.cbsub;
e = sqrt(diag(cov))';

subplot(1,4,1)
loglog(rsub_arr, d, 'k.', 'markersize', 10);hold on
loglog(rsub_arr, e,'b.-');
loglog(rsub_arr, sqrt(1./diag(covi)),'r.-');
h=legend({'Excess','sqrt(diag(Excess Cov))','1/sqrt(diag(Excess Cov^{-1}))'},...
        'Location','northeast');
set(h,'fontsize',10)
loglog(rsub_arr, -d, 'ko', 'markersize', 5);hold on
xlim([4e-1,1.1e3])
ylim([1e-1,1e2])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('I [nW/m^2/sr]', 'fontsize',15);


diff_list = [];

diff = e;
diff_list = [diff_list;diff];

diff = e;
diff(2:2:end) = -diff(2:2:end);
diff_list = [diff_list;diff];

diff = e.*(0.1:0.1:1.5);
diff_list = [diff_list;diff];

diff = e.*(0.1:0.1:1.5);
diff(2:2:end) = -diff(2:2:end);
diff_list = [diff_list;diff];

diff = e.*(1.5:-0.1:0.1);
diff_list = [diff_list;diff];

diff = e.*(1.5:-0.1:0.1);
diff(2:2:end) = -diff(2:2:end);
diff_list = [diff_list;diff];

subplot(1,4,2)
for i=1:2
    diff = diff_list(i,:);
    chi2indept = sum(diff.^2./e.^2);
    chi2full = (diff*covi*diff');
    chi2mat = (diff'*diff).*covi;
    semilogx(rsub_arr, diff./e,'-','DisplayName',...
        strcat('\chi^2_{indept} =',num2str(chi2indept,'%.2f'),...
        ', \chi^2_{full} =',num2str(chi2full,'%.2f')));hold on
end

h=legend('show','Location','southeast');
set(h,'fontsize',10)
xlim([4e-1,1.1e3])
ylim([-2,2])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('\Delta/\sigma', 'fontsize',15);

subplot(1,4,3)
for i=3:4
    diff = diff_list(i,:);
    chi2indept = sum(diff.^2./e.^2);
    chi2full = (diff*covi*diff');
    semilogx(rsub_arr, diff./e,'-','DisplayName',...
        strcat('\chi^2_{indept} =',num2str(chi2indept,'%.2f'),...
        ', \chi^2_{full} =',num2str(chi2full,'%.2f')));hold on
end
h=legend('show','Location','southeast');
set(h,'fontsize',10)
xlim([4e-1,1.1e3])
ylim([-2,2])
ylabel('\Delta/\sigma', 'fontsize',15);

subplot(1,4,4)
for i=5:6
    diff = diff_list(i,:);
    chi2indept = sum(diff.^2./e.^2);
    chi2full = (diff*covi*diff');
    semilogx(rsub_arr, diff./e,'-','DisplayName',...
        strcat('\chi^2_{indept} =',num2str(chi2indept,'%.2f'),...
        ', \chi^2_{full} =',num2str(chi2full,'%.2f')));hold on
end
h=legend('show','Location','southeast');
set(h,'fontsize',10)
xlim([4e-1,1.1e3])
ylim([-2,2])
ylabel('\Delta/\sigma', 'fontsize',15);
%% chi2 matrix

ifield = 4;
masklim = false;
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackdat_%s',...
        loaddir,dt.name),'stackdatall');

figure
setwinsize(gcf,1500,600)
im = 2;

stackdat = stackdatall(im).stackdat;
m_min = stackdat.m_min;
m_max = stackdat.m_max;
rsub_arr = stackdat.rsub_arr;

cov = stackdat.cov.cbsub;
covi = stackdat.cov.cbsub_inv;
d = stackdat.excess.cbsub;
e = sqrt(diag(cov))';

diff_list = [];
diff = e;
diff_list = [diff_list;diff];
diff = e;
diff(2:2:end) = -diff(2:2:end);
diff_list = [diff_list;diff];

for i=1:2
    diff = diff_list(i,:);
    chi2indept = sum(diff.^2./e.^2);
    chi2full = (diff*covi*diff');
    chi2mat = (diff'*diff).*covi;
    subplot(2,3,1)
    semilogx(rsub_arr, diff./e,'-');hold on
    
    subplot(2,3,2)
    semilogx(rsub_arr, diff,'-','DisplayName',strcat('case',num2str(i),...
        ': \chi^2_{indept} =',num2str(chi2indept,'%.2f'),...
        ', \chi^2_{full} =',num2str(chi2full,'%.2f')));hold on 
    
    subplot(2,3,3+i)
    imageclip(chi2mat);
    caxis([-7,7]);
    title(strcat('(Case ',num2str(i),') \Delta_i C^{-1}_{ij} \Delta_j'));
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)

    subplot(2,3,6)
    cum_chi2 = zeros(size(rsub_arr));
    for irsub=1:numel(rsub_arr)
        cum_chi2(irsub) = sum(sum(chi2mat(1:irsub,1:irsub)))/irsub;
    end
    loglog(rsub_arr, cum_chi2,'.-','DisplayName',strcat('case',num2str(i)));hold on
    
end

subplot(2,3,3)
imageclip(covi);
caxis([-1.5,1.5]);
h = colorbar;
ylabel(h, '[nW/m^2/sr]^{-2}');
title('C^{-1}');
xticks([3:4:15])
xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
    num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
xtickangle(45)
yticks([3:4:15])
yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
    num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
ytickangle(45)
if im==1
    ylabel('Excess Cov Inv.', 'fontsize',15);
end

subplot(2,3,1)
xlim([4e-1,1.1e3])
ylim([-1.2,1.2])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('\Delta/\sigma', 'fontsize',15);

subplot(2,3,2)
h=legend('show','Location','northeast');
set(h,'fontsize',10)
xlim([4e-1,1.1e3])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('\Delta [nW/m^2/sr]', 'fontsize',12);

subplot(2,3,6)
h=legend('show','Location','northeast');
set(h,'fontsize',10)
xlim([4e-1,1.1e3])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('\Sigma \Delta C^{-1} \Delta (<r) / DoF (<r)', 'fontsize',12);
%% plot fit with fiducial value
ifield = 4;
masklim = false;
dx = 1200;
nbins = 25;
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackdat_%s',...
        loaddir,dt.name),'stackdatall');
   
figure
setwinsize(gcf,1200,500)
for im = 1:4
m_min = stackdatall(im).stackdat.m_min;
m_max = stackdatall(im).stackdat.m_max;
rsub_arr = stackdatall(im).stackdat.rsub_arr;
r_arr = stackdatall(im).stackdat.r_arr;

% data excess
dat = stackdatall(im).stackdat.all.profcbg;
edat = stackdatall(im).stackdat.excess.cbsub;
cov_inv = stackdatall(im).stackdat.cov.cbsub_inv;
cov_mat = stackdatall(im).stackdat.cov.cbsub;

% hsc clustering excess
weight = stackdatall(im).stackdat.all.profhitg;
hscclusdat = get_hsc_clus_prof(flight, inst, masklim, weight);
clus = hscclusdat(4).linefitmodelsub;

% PSF x window map
psfwin_arr = stackdatall(im).stackdat.bgsub.profcbpsf;
radmap = make_radius_map(zeros(2*dx+1),dx+1,dx+1).*0.7;
psfwin_map = spline(r_arr(1:18),psfwin_arr(1:18),radmap);
psfwin_map(radmap > r_arr(18)) = 0;
profile = radial_prof(psfwin_map,ones(2*dx+1),dx+1,dx+1,1,nbins);
psfwin_arr = profile.prof./profile.prof(1);
[psfwinsub_arr] = profile_radial_binning(psfwin_arr,weight,1);

% model map
rrmap = make_radius_map(zeros(201),101,101).*0.7;
[hsc_map,hsc_map_sub,R200,param_arr] = HSC_Wang19_prof(rrmap,im,true);

hsc_map1 = squeeze(hsc_map_sub(1,:,:));
psfwinhsc_map1 = conv2(psfwin_map,hsc_map1,'same');

hsc_map2 = squeeze(hsc_map_sub(2,:,:));
psfwinhsc_map2 = conv2(psfwin_map,hsc_map2,'same');

Ie1  = param_arr(1,1);
n1 = param_arr(1,2);
xe1 = param_arr(1,3);
Ie2  = param_arr(2,1);
n2 = param_arr(2,2);
xe2 = param_arr(2,3);
chi2best = 1e8;
subplot(2,4,im)
hsc_map2i = sersic(rrmap./R200,Ie2,n2,xe1+xe2);
psfwinhsc_map2i = conv2(psfwin_map, hsc_map2i,'same');
profile = radial_prof(psfwinhsc_map2i + psfwinhsc_map1,...
    ones(2*dx+1),dx+1,dx+1,1,nbins);
prof = profile.prof / profile.prof(1);
efull = prof.*dat(1) - psfwin_arr.*dat(1);
esub = profile_radial_binning(efull,weight,1);
clus_scale = 1;
m = clus_scale.*clus + esub;
diff = edat - m;
chi2 = diff*cov_inv*diff';
chi2mat = (diff'*diff).*cov_inv;
fprintf('xe = %.3f, clus_scale = %d,chi2 = %.2e\n',...
    xe2,clus_scale,chi2); 
loglog(rsub_arr,m,'DisplayName',...
sprintf('(xe2,A_{clus}) = (%.4f, %d) \n chi^2 = %.1f',...
xe2,clus_scale,chi2));hold on

h=legend('show','Location','southwest');
set(h,'fontsize',10)
legend boxoff
loglog(rsub_arr,edat,'k.','markersize',10);hold on
loglog(rsub_arr,-edat,'ko','markersize',5);
errorbar(rsub_arr,edat,sqrt(diag(cov_mat))','k.','markersize',10);
errorbar(rsub_arr,-edat,sqrt(diag(cov_mat))','ko','markersize',5);
xlim([4e-1,1.1e3])
ylim([1e-2,1e2])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('Excess I [nW/m^2/sr]', 'fontsize',12);
title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);

subplot(2,4,im+4)
imageclip(chi2mat);
xticks([3:4:15])
xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
    num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
xtickangle(45)
yticks([3:4:15])
yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
    num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
ytickangle(45)
title(strcat(' \Delta_i C^{-1}_{ij} \Delta_j'));
end
%%
ifield = 4;
masklim = false;
dx = 1200;
nbins = 25;
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackdat_%s',...
        loaddir,dt.name),'stackdatall');
   
figure
setwinsize(gcf,1200,500)

for im = 1:4
m_min = stackdatall(im).stackdat.m_min;
m_max = stackdatall(im).stackdat.m_max;
rsub_arr = stackdatall(im).stackdat.rsub_arr;
r_arr = stackdatall(im).stackdat.r_arr;

% data excess
dat = stackdatall(im).stackdat.all.profcbg;
edat = stackdatall(im).stackdat.excess.cbsub;
cov_inv = stackdatall(im).stackdat.cov.cbsub_inv;
cov_mat = stackdatall(im).stackdat.cov.cbsub;

% hsc clustering excess
weight = stackdatall(im).stackdat.all.profhitg;
hscclusdat = get_hsc_clus_prof(flight, inst, masklim, weight);
clus = hscclusdat(4).linefitmodelsub;

% PSF x window map
psfwin_arr = stackdatall(im).stackdat.bgsub.profcbpsf;
radmap = make_radius_map(zeros(2*dx+1),dx+1,dx+1).*0.7;
psfwin_map = spline(r_arr(1:18),psfwin_arr(1:18),radmap);
psfwin_map(radmap > r_arr(18)) = 0;
profile = radial_prof(psfwin_map,ones(2*dx+1),dx+1,dx+1,1,nbins);
psfwin_arr = profile.prof./profile.prof(1);
[psfwinsub_arr] = profile_radial_binning(psfwin_arr,weight,1);

% model map
rrmap = make_radius_map(zeros(201),101,101).*0.7;
[hsc_map,hsc_map_sub,R200,param_arr] = HSC_Wang19_prof(rrmap,im,true);

hsc_map1 = squeeze(hsc_map_sub(1,:,:));
psfwinhsc_map1 = conv2(psfwin_map,hsc_map1,'same');

hsc_map2 = squeeze(hsc_map_sub(2,:,:));
psfwinhsc_map2 = conv2(psfwin_map,hsc_map2,'same');

Ie1  = param_arr(1,1);
n1 = param_arr(1,2);
xe1 = param_arr(1,3);
Ie2  = param_arr(2,1);
n2 = param_arr(2,2);
xe2 = param_arr(2,3);
xe2_arr = [0.01,0.02,0.03,0.04,0.05];
clus_scale_arr = [1,2,3];
chi2_arr = zeros([numel(xe2_arr,clus_scale_arr)]);
subplot(2,4,im)
for ixe2= 1:numel(xe2_arr)
    xe2 = xe2_arr(ixe2);
    hsc_map2i = sersic(rrmap./R200,Ie2,n2,xe1+xe2);
    psfwinhsc_map2i = conv2(psfwin_map, hsc_map2i,'same');
    profile = radial_prof(psfwinhsc_map2i + psfwinhsc_map1,...
        ones(2*dx+1),dx+1,dx+1,1,nbins);
    prof = profile.prof / profile.prof(1);
    efull = prof.*dat(1) - psfwin_arr.*dat(1);
    esub = profile_radial_binning(efull,weight,1);
    for iclus_scale = 1:numel(clus_scale_arr)
        clus_scale = clus_scale_arr(iclus_scale);
        m = clus_scale.*clus + esub;
        diff = edat - m;
        chi2 = diff*cov_inv*diff';
        chi2_arr(ixe2,iclus_scale) = chi2;
        fprintf('xe = %.3f, clus_scale = %d,chi2 = %.2e\n',...
            xe2,clus_scale,chi2);
        if xe2==0.01
            loglog(rsub_arr,clus_scale.*clus,'DisplayName',...
            sprintf('A_{clus} = %d',clus_scale));hold on
        end        
    end
    loglog(rsub_arr,esub,'DisplayName',...
    sprintf('xe2 = %.2f',xe2));hold on
end
if im==4
    h=legend('show','Location','northeast');
    set(h,'fontsize',10)
    legend boxoff
end
loglog(rsub_arr,edat,'k.','markersize',10);hold on
loglog(rsub_arr,-edat,'ko','markersize',5);
errorbar(rsub_arr,edat,sqrt(diag(cov_mat))','k.','markersize',10);
errorbar(rsub_arr,-edat,sqrt(diag(cov_mat))','ko','markersize',5);
xlim([4e-1,1.1e3])
ylim([1e-2,1e2])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('Excess I [nW/m^2/sr]', 'fontsize',12);
title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);

subplot(2,4,im+4)
for iclus_scale=1:numel(clus_scale_arr)
    clus_scale = clus_scale_arr(iclus_scale);
    plot(xe2_arr, chi2_arr(:,iclus_scale),'o-','DisplayName',...
            sprintf('A_{clus} = %d',clus_scale));hold on
end
xlabel('xe2', 'fontsize',15);
ylabel('\chi^2', 'fontsize',12);
title(strcat('\chi^2'));
if im==4
    h=legend('show','Location','northwest');
    set(h,'fontsize',10)
    legend boxoff
end

end
