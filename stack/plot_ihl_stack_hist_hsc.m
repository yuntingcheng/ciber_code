function plot_ihl_stack_hist_hsc(flight,inst,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addOptional('f_ihl',0,@isnumeric);
  p.addOptional('rvir',1,@isnumeric);
  p.addOptional('rmin',nan,@isnumeric);
  
  p.parse(flight,inst,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  f_ihl    = p.Results.f_ihl;
  rvir     = p.Results.rvir;
  rmin     = p.Results.rmin;
  
  clear p varargin;

%%
mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
%%
for hsc_idx=[0,1,2,3]
    if hsc_idx==0
        field='SIDES';
    else
        field='HSC';
    end
    
    if strcmp(field,'SIDES')
        name = 'SIDES';
        if rvir==1
            load(sprintf('%s/ihlprofdatsim_hist%d',loaddir,f_ihl*100),...
                'ihlprofdat');
        else
            load(sprintf('%s/ihlprofdatsim_hist%d_rv%d',...
                loaddir,f_ihl*100,rvir),'ihlprofdat');
        end
    elseif strcmp(field,'HSC')
        name = HSC_fields_info(hsc_idx);
        if rvir==1
            load(sprintf('%s/ihlprofdathsc_%s_hist%d',loaddir,name,f_ihl*100),...
                'ihlprofdat');
        else
            load(sprintf('%s/ihlprofdathsc_%s_hist%d_rv%d',...
                name,loaddir,f_ihl*100,rvir),'ihlprofdat');
        end
    end

figure
setwinsize(gcf,1200,1200)
for im=1:3
    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;
    counts = ihlprofdat.data(im).counts;
    countg = ihlprofdat.data(im).countg;
    if numel(counts)==0
        continue
    end
    r_arr = ihlprofdat.r_arr;

    profscb = ihlprofdat.data(im).profscb;
    profscb_err = ihlprofdat.data(im).profscb_err;
    profsps = ihlprofdat.data(im).profsps;
    profsps_err = ihlprofdat.data(im).profsps_err;
    profgcb = ihlprofdat.data(im).profgcb;
    profgcb_err = ihlprofdat.data(im).profgcb_err;
    profgps = ihlprofdat.data(im).profgps;
    profgps_err = ihlprofdat.data(im).profgps_err;

    profscbbk = ihlprofdat.bk(im).profscb;
    profscb_errbk = ihlprofdat.bk(im).profscb_err;
    profspsbk = ihlprofdat.bk(im).profsps;
    profsps_errbk = ihlprofdat.bk(im).profsps_err;

    subplot(3,3,im)
    semilogx(r_arr.*0.98,profscb,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgcb,'b','markersize',10);
    plot(r_arr.*0.98,profsps,'m','markersize',10);
    plot(r_arr.*1.02,profgps,'c','markersize',10);
    h=legend({'CIBER stars','CIBER galaxies','Sim stars','Sim galaxies'},...
        'Location','northeast');
    set(h,'fontsize',10)
    legend boxoff
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    errorbar(r_arr.*0.98,profsps,profsps_err,'m.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'c.','markersize',10);
    xlim([4e-1,50])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
      ' (',num2str(counts),' stars, ',num2str(countg), ' galaxies)'),'fontsize',15)
            
    subplot(3,3,im + 3)
    semilogx(r_arr.*0.98,profscb,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgcb,'b','markersize',10);
    plot(r_arr,profscbbk,'k','markersize',10);
    h=legend({'CIBER stars','CIBER galaxies','background'},...
        'Location','northwest');
    set(h,'fontsize',10)
    legend boxoff
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    errorbar(r_arr,profscbbk,profscb_errbk,'k.','markersize',10);
    xlim([4e-1,1e3])
    ylim([-1,10])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)

    subplot(3,3,im + 6)
    semilogx(r_arr.*0.98,profsps,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgps,'b','markersize',10);
    plot(r_arr,profspsbk,'k','markersize',10);
    h=legend({'Sim stars','Sim galaxies'},...
        'Location','northwest');
    set(h,'fontsize',10)
    legend boxoff
    errorbar(r_arr.*0.98,profsps,profsps_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'b.','markersize',10);
    errorbar(r_arr,profspsbk,profsps_errbk,'k.','markersize',10);
    xlim([4e-1,1e3])
    ylim([-0.05,0.1])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end
s=suptitle(strrep(name,'_','\_'));
set(s,'FontSize',15);
savename=strcat(pltsavedir,name,'_stackprof_hist_sim',num2str(f_ihl*100));
% print(savename,'-dpng');%close
end
%% plot the excess profile%%%
for hsc_idx=[0,1,2,3]
    if hsc_idx==0
        field='SIDES';
    else
        field='HSC';
    end
    
    if strcmp(field,'SIDES')
        name = 'SIDES';
        if rvir==1
            load(sprintf('%s/ihlprofdatsim_hist%d',loaddir,f_ihl*100),...
                'ihlprofdat');
        else
            load(sprintf('%s/ihlprofdatsim_hist%d_rv%d',...
                loaddir,f_ihl*100,rvir),'ihlprofdat');
        end
    elseif strcmp(field,'HSC')
        name = HSC_fields_info(hsc_idx);
        if rvir==1
            load(sprintf('%s/ihlprofdathsc_%s_hist%d',loaddir,name,f_ihl*100),...
                'ihlprofdat');
        else
            load(sprintf('%s/ihlprofdathsc_%s_hist%d_rv%d',...
                name,loaddir,f_ihl*100,rvir),'ihlprofdat');
        end
    end

figure
setwinsize(gcf,1200,360)
for im=1:3
    m_min = ihlprofdat.data(im).m_min;
    m_max = ihlprofdat.data(im).m_max;
    counts = ihlprofdat.data(im).counts;
    countg = ihlprofdat.data(im).countg;
    if numel(counts)==0
        continue
    end

    r_arr = ihlprofdat.r_arr;

    ecb = ihlprofdat.excess(im).diffcb;
    ecb_err = ihlprofdat.excess(im).diffcb_err1;
    eps = ihlprofdat.excess(im).diffps;
    eps_err = ihlprofdat.excess(im).diffps_err1;
    ihl = ihlprofdat.excess(im).diff;
    ihl_err = ihlprofdat.excess(im).diff_err1;
    
    subplot(1,3,im)
    loglog(r_arr.*0.98,ecb,'r.','markersize',10);hold on
    loglog(r_arr.*1.02,eps,'b.','markersize',10);hold on
    loglog(r_arr,ihl,'k.','markersize',10);hold on
    
    sp = find(r_arr > 100);
    v = ihlprofdat.excess(im).diff100;
    e = ihlprofdat.excess(im).diff_err100;
    if v>=0 & v-e >=0
        fill([r_arr(sp),flip(r_arr(sp))],...
            [(v+e) * ones(size(sp)),(v-e) * ones(size(sp))],...
                [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
        plot(r_arr(sp),v*ones(size(sp)),'color',[0.2,0.4,0.2],'linewidth',2);
    elseif v>=0 & v-e<0
        fill([r_arr(sp),flip(r_arr(sp))],...
            [(v+e) * ones(size(sp)),1e-8 * ones(size(sp))],...
                [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
        plot(r_arr(sp),v*ones(size(sp)),'color',[0.2,0.4,0.2],'linewidth',2);
    elseif v<0 & v+e>0
        fill([r_arr(sp),flip(r_arr(sp))],...
            [(v+e) * ones(size(sp)),1e-8 * ones(size(sp))],...
                [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
        plot(r_arr(sp),-v*ones(size(sp)),'--','color',[0.2,0.4,0.2],'linewidth',2);
    elseif v<0 & v+e<=0
        fill([r_arr(sp),flip(r_arr(sp))],...
            [-(v-e) * ones(size(sp)),-(v+e) * ones(size(sp))],...
                [0.6,0.4,0.2],'facealpha',0.5,'EdgeColor','none');
        plot(r_arr(sp),-v*ones(size(sp)),'--','color',[0.6,0.4,0.2],'linewidth',2);
    end

    h=legend({'Excess all(bright+faint+IHL)','Excess bright', ...
        'Excess faint+IHL','Excess faint+IHL> 100 arcsec'},'Location','northeast');
    set(h,'fontsize',10)
    legend boxoff

    eps_err_low = eps_err;
    for ii=1:numel(eps)
        if eps(ii)>0 & eps_err(ii) > eps(ii)
            eps_err_low(ii) = eps(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*1.02,eps,eps_err_low,eps_err,'b.','markersize',10);
    loglog(r_arr.*1.02,-eps,'bo','markersize',5);hold on
    eps_err_low = eps_err;
    for ii=1:numel(eps)
        if eps(ii)<0 & eps_err(ii) > abs(eps(ii))
            eps_err_low(ii) = -eps(ii) - 1e-8;
        end
    end    
    errorbar(r_arr.*1.02,-eps,eps_err_low,eps_err,'bo','markersize',5);
    
    ecb_err_low = ecb_err;
    for ii=1:numel(eps)
        if ecb(ii)>0 & ecb_err(ii) > ecb(ii)
            ecb_err_low(ii) = ecb(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*0.98,ecb,ecb_err_low,ecb_err,'r.','markersize',10);
    loglog(r_arr.*0.98,-ecb,'ro','markersize',5);
    ecb_err_low = ecb_err;
    for ii=1:numel(ecb)
        if ecb(ii)<0 & ecb_err(ii) > abs(ecb(ii))
            ecb_err_low(ii) = -ecb(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*0.98,-ecb,ecb_err_low,ecb_err,'ro','markersize',5);
    
    ihl_err_low = ihl_err;
    for ii=1:numel(ihl)
        if ihl(ii)>0 & ihl_err(ii) > ihl(ii)
            ihl_err_low(ii) = ihl(ii) - 1e-8;
        end
    end
    errorbar(r_arr,ihl,ihl_err_low,ihl_err,'k.','markersize',10);
    loglog(r_arr,-ihl,'ko','markersize',5);
    
    ihl_err_low = ihl_err;
    for ii=1:numel(ecb)
        if ihl(ii)<0 & ihl_err(ii) > abs(ihl(ii))
            ihl_err_low(ii) = -ihl(ii) - 1e-8;
        end
    end    
    errorbar(r_arr,-ihl,ihl_err_low,ihl_err,'ko','markersize',5);
    xlim([4e-1,1e3])
    ylim([1e-4,1e4])
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
      ' (',num2str(counts),' stars, ',num2str(countg), ' galaxies)'),'fontsize',15)
    xlabel('arcsec', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)

end

suptitle(strrep(name,'_','\_'));
savename=strcat(pltsavedir,name,'_excessprof_hist');
print(savename,'-dpng');%close
end

%% plot the excess profile%%%

% get the average access across field
for im=1:3
    ihldat = zeros([5,25]);
    errdat = zeros([5,25]);
    for ifield = 4:8
        dt=get_dark_times(flight,inst,ifield);
        loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
        if rmin==2
            load(sprintf('%s/%s_ihlprofdat_hist_rmin2',loaddir,dt.name),...
                'ihlprofdat');
        elseif isnan(rmin)
            load(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');
        end
        ihl = ihlprofdat.excess(im).diff;
        ihl_err = ihlprofdat.excess(im).diff_err1;
        ihldat(ifield-3,:) = ihl;
        errdat(ifield-3,:) = ihl_err;
    end
    ihltot = sum(ihldat./errdat.^2)./sum(1./errdat.^2);
    ihltot_err = sqrt(1./(sum(1./errdat.^2)));
    ihlcbdat(im).ihltot = ihltot;
    ihlcbdat(im).ihltot_err = ihltot_err;
end
r_arr = ihlprofdat.r_arr;

figure
setwinsize(gcf,1200,900)
c_arr = {'k','b','r','m','c'};
off_arr=[1,0.98,0.96,1.02,1.04];
count=0;
legendname = {};
ihlavg = zeros(3,numel(ihltot));
ihl_erravg = zeros(3,numel(ihltot));
for hsc_idx=[-999,0,1,2,3]
    count=count+1;
    c = c_arr{count};
    off = off_arr(count);

    if hsc_idx==-999
        field='CIBER';
        legendname{end+1} = field;
    elseif hsc_idx==0
        field='SIDES';
        legendname{end+1} = field;
    else
        field='HSC';
    end
    
    if strcmp(field,'SIDES')
        name = 'SIDES';
        if rvir==1
            load(sprintf('%s/ihlprofdatsim_hist%d',loaddir,f_ihl*100),...
                'ihlprofdat');
        else
            load(sprintf('%s/ihlprofdatsim_hist%d_rv%d',...
                loaddir,f_ihl*100,rvir),'ihlprofdat');
        end
    elseif strcmp(field,'HSC')
        name = HSC_fields_info(hsc_idx);
        legendname{end+1} = strcat(field,'-',strrep(name,'_','\_'));
        if rvir==1
            load(sprintf('%s/ihlprofdathsc_%s_hist%d',loaddir,name,f_ihl*100),...
                'ihlprofdat');
        else
            load(sprintf('%s/ihlprofdathsc_%s_hist%d_rv%d',...
                name,loaddir,f_ihl*100,rvir),'ihlprofdat');
        end
    end

    for im=1:3
        m_min = ihlprofdat.data(im).m_min;
        m_max = ihlprofdat.data(im).m_max;
        if hsc_idx~=-999
            counts = ihlprofdat.data(im).counts;
            if numel(counts)==0
                continue
            end
            ihl = ihlprofdat.excess(im).diff;
            ihl_err = ihlprofdat.excess(im).diff_err1;
        else
            ihl = ihlcbdat(im).ihltot;
            ihl_err = ihlcbdat(im).ihltot_err;
        end

        subplot(3,3,im+3)
        semilogx(r_arr*off,ihl./ihlcbdat(im).ihltot,...
            '.','color',c,'markersize',10);hold on
        errorbar(r_arr*off,ihl./ihlcbdat(im).ihltot,...
            ihl_err./ihlcbdat(im).ihltot,'.','color',c,'markersize',10);
        xlim([4e-1,1e3])
        ylim([-0.5,3])
        xlabel('arcsec', 'fontsize',15)
        ylabel('I/I_{CIBER}', 'fontsize',15)

        if hsc_idx==-999
            subplot(3,3,im+6)
            semilogx(r_arr*off,ihl./ihlcbdat(im).ihltot,...
                '.','color',c,'markersize',10);hold on
            errorbar(r_arr*off,ihl./ihlcbdat(im).ihltot,...
                ihl_err./ihlcbdat(im).ihltot,'.','color',c,'markersize',10);
        else
            ihlavg(im,:) = ihlavg(im,:) + (ihl./ihl_err.^2);
            ihl_erravg(im,:) = ihl_erravg(im,:) + (1./ihl_err.^2);
        end
        
        subplot(3,3,im)
        hh(count) = loglog(r_arr*off,ihl,'.','color',c,'markersize',10);hold on

        ihl_err_low = ihl_err;
        for ii=1:numel(ihl)
            if ihl(ii)>0 & ihl_err(ii) > ihl(ii)
                ihl_err_low(ii) = ihl(ii) - 1e-8;
            end
        end
        errorbar(r_arr*off,ihl,ihl_err_low,ihl_err,'.','color',c,'markersize',10);
        loglog(r_arr*off,-ihl,'o','color',c,'markersize',5);

        ihl_err_low = ihl_err;
        for ii=1:numel(ihl)
            if ihl(ii)<0 & ihl_err(ii) > abs(ihl(ii))
                ihl_err_low(ii) = -ihl(ii) - 1e-8;
            end
        end    
        errorbar(r_arr*off,-ihl,ihl_err_low,ihl_err,'o','color',c,'markersize',5);
        xlim([4e-1,1e3])
        ylim([1e-4,1e2])
        title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15)
        xlabel('arcsec', 'fontsize',15)
        ylabel('I [nW/m^2/sr]', 'fontsize',15)
    end
end   
h=legend(hh,legendname);
    set(h,'fontsize',7,'Location','northeast')
    legend boxoff

for im=1:3
    ihl = ihlavg(im,:)./ihl_erravg(im,:);
    ihl_err = sqrt(1./ihl_erravg(im,:));
    subplot(3,3,im+6)
    hh=semilogx(r_arr*off,ihl./ihlcbdat(im).ihltot,...
        'b.','markersize',10);hold on
    errorbar(r_arr*off,ihl./ihlcbdat(im).ihltot,...
        ihl_err./ihlcbdat(im).ihltot,'b.','markersize',10);
    xlim([4e-1,1e3])
    ylim([-0.5,3])
    xlabel('arcsec', 'fontsize',15)
    ylabel('I/I_{CIBER}', 'fontsize',15)
    h=legend(hh,{'sims average'});
    set(h,'fontsize',10,'Location','northeast')
    legend boxoff
end
        
savename=strcat(pltsavedir,'HSC_excessprof_hist');
print(savename,'-dpng');%close
%%
return