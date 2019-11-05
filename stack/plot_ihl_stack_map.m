function plot_ihl_stack_map(flight,inst,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addOptional('masklim',false,@islogical);
  p.addOptional('savefig',false,@islogical);
  
  p.parse(flight,inst,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  masklim = p.Results.masklim;
  savefig = p.Results.savefig;
  
  clear p varargin;

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
    errl = err;
    errl(errl > abs(prof)) = abs(prof(errl > abs(prof))) - 1e-10;
    loglog(r_arr, prof,'k.','markersize',10);hold on
    errorbar(r_arr, prof,errl,err,'k.','markersize',10);
    errorbar(r_arr, -prof,errl,err,'ko','markersize',5);
    prof = excessdat.excess.profpsg;
    err = sqrt(diag(excessdat.excov.covps))';
    errl = err;
    errl(errl > abs(prof)) = abs(prof(errl > abs(prof))) - 1e-10;
    errorbar(r_arr.*1.03, prof,errl,err,'b.','markersize',10);
    errorbar(r_arr.*1.03, -prof,errl,err,'bo','markersize',5);
    prof = hscclusdat(im).prof;
    err = hscclusdat(im).err;
    errl = err;
    errl(errl > abs(prof)) = abs(prof(errl > abs(prof))) - 1e-10;
    errorbar(r_arr.*0.97, prof,errl,err,'r.','markersize',10);
    errorbar(r_arr.*0.97, -prof,errl,err,'ro','markersize',5);
    
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

for itype=1:4
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
    case 4
        exdat = excessdat.excov;
        name = 'excess';
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
    caxis([-1,1])
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

%% plot the diag of cov
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
setwinsize(gcf,1400,250)
for im=1:4
    stackdat = stackdatall(im).stackdat;
    excessdat = excessdatall(im).excessdat;
    m_min = stackdat.m_min;
    m_max = stackdat.m_max;
    r_arr = stackdat.r_arr;

    subplot(1,4,im)
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
        h=legend({'data err','BG err','PSF err', 'Excess err',...
            'Excess err (jackknife)', 'EExcess err (jackknife w/ norm)',...
            'Excess Data'},'Location','northeast');
        set(h,'fontsize',6)
    end
    
    loglog(r_arr, -dat, 'ko', 'markersize', 5);hold on
    xlim([4e-1,1.1e3])
    ylim([1e-1,1.1e2])
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return