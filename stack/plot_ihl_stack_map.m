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

%%
mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

%% compile excess from all fields

for ifield = 4:8
    dt=get_dark_times(flight,inst,ifield);

    loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
    if masklim
        load(sprintf('%s/stackdat_%s_masklim',...
                loaddir,dt.name),'stackdatall');        
    else
        load(sprintf('%s/stackdat_%s',...
                loaddir,dt.name),'stackdatall');
    end
    for im=1:4
        stackdat = stackdatall(im).stackdat;
        excess(ifield).r_arr = stackdat.r_arr;
        excess(ifield).mag(im).m_min = stackdat.m_min;
        excess(ifield).mag(im).m_max = stackdat.m_max;
        excess(ifield).count(im).counts = stackdat.all.counts;
        excess(ifield).count(im).countg = stackdat.all.countg;
        excess(ifield).cb(im).ihl = stackdat.excess.diffcbpsf;
        excess(ifield).cb(im).err = stackdat.excess.diffcbpsf_err;
        excess(ifield).cb(im).ihl100 = stackdat.excess.diffcbpsf100;
        excess(ifield).cb(im).err100 = stackdat.excess.diffcbpsf_err100;
        excess(ifield).ps(im).ihl = stackdat.excess.diffpspsf;
        excess(ifield).ps(im).err = stackdat.excess.diffpspsf_err;

        % src profile
        profsrc = stackdatall(2).stackdat.psfcb;
        profsrc_err = stackdatall(2).stackdat.psfcb_err;
        norm = stackdat.norm.profcbg(1)/profsrc(1);
        excess(ifield).cb(im).src = profsrc.*norm;
        excess(ifield).cb(im).src_err = profsrc_err.*norm;
        
        % fluc error
        errfluc = [];
        errfluc100 = [];
        for iter=1:50
            try
                if masklim
                    load(sprintf('%s/fluc/stackdatfluc_%s_masklim_iter%d',...
                            loaddir,dt.name,iter),'stackdatfluc');        
                else
                    load(sprintf('%s/fluc/stackdatfluc_%s_iter%d',...
                            loaddir,dt.name,iter),'stackdatfluc');
                end
                errfluc = [errfluc;...
                    stackdatfluc(im).stackdat.excess.diff];
                errfluc100 = [errfluc100;...
                    stackdatfluc(im).stackdat.excess.diff100];
            catch
                continue
            end
        end
        
        if numel(errfluc100)<=1
            excess(ifield).cb(im).errf = 0;
        else
            excess(ifield).cb(im).errf = nanstd(errfluc);
        end
        
        if numel(errfluc100)<=1
            excess(ifield).cb(im).errf100 = 0;
        else
            excess(ifield).cb(im).errf100 = nanstd(errfluc100);
        end
        excess(ifield).cb(im).errdf = sqrt(excess(ifield).cb(im).err.^2 + ...
            excess(ifield).cb(im).errf.^2);
        excess(ifield).cb(im).errdf100 = sqrt(excess(ifield).cb(im).err100.^2 + ...
            excess(ifield).cb(im).errf100.^2);        
    end
end

% field avg %
r_arr = stackdat.r_arr;
for im=1:4
ihlcbdat = zeros([5,numel(r_arr)]);
errcbdat = zeros([5,numel(r_arr)]);
ihlpsdat = zeros([5,numel(r_arr)]);
errpsdat = zeros([5,numel(r_arr)]);
srccbdat = zeros([5,numel(r_arr)]);
srccbdat_err = zeros([5,numel(r_arr)]);

cbihldat100 = zeros([5,1]);
cberrdat100 = zeros([5,1]);
for ifield = 4:8
    ihlcbdat(ifield-3,:) = excess(ifield).cb(im).ihl;
    errcbdat(ifield-3,:) = excess(ifield).cb(im).err;
    ihlpsdat(ifield-3,:) = excess(ifield).ps(im).ihl;
    errpsdat(ifield-3,:) = excess(ifield).ps(im).err;
    srccbdat(ifield-3,:) = excess(ifield).cb(im).src;
    srccbdat_err(ifield-3,:) = excess(ifield).cb(im).src_err;
    cbihldat100(ifield-3,:) = excess(ifield).cb(im).ihl100;
    cberrdat100(ifield-3,:) = excess(ifield).cb(im).err100;

end
ihltotcb = sum(ihlcbdat./errcbdat.^2)./sum(1./errcbdat.^2);
ihltotcb_err = sqrt(1./(sum(1./errcbdat.^2)));
excess(1).cb(im).ihl = ihltotcb;
excess(1).cb(im).err = ihltotcb_err;
ihltotps = sum(ihlpsdat./errpsdat.^2)./sum(1./errpsdat.^2);
ihltotps_err = sqrt(1./(sum(1./errpsdat.^2)));
excess(1).ps(im).ihl = ihltotps;
excess(1).ps(im).err = ihltotps_err;

srctotcb = sum(srccbdat./srccbdat_err.^2)./sum(1./srccbdat_err.^2);
srctotcb_err = sqrt(1./(sum(1./srccbdat_err.^2)));
excess(1).cb(im).src = srctotcb;
excess(1).cb(im).src_err = srctotcb_err;

excess(1).cb(im).ihl100 = sum(cbihldat100./cberrdat100.^2)./sum(1./cberrdat100.^2);
excess(1).cb(im).err100 = sqrt(1./sum(1./cberrdat100.^2));

end

% hsc
[psf_arr,~,~] =  PSF_stacked_profile(flight,inst,8);
for im=1:4
    excb_all = zeros([12,numel(r_arr)]);
    exps_all = zeros([12,numel(r_arr)]);
    excb_err_all = zeros([12,numel(r_arr)]);
    exps_err_all = zeros([12,numel(r_arr)]);
    for hsc_idx=0:11
        [name,~] = HSC_fields_info(hsc_idx);
        loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
        if masklim
            load(sprintf('%s/hsc/stackdathsc_%s_masklim',...
                loaddir,name),'stackdathsc');  
            load(sprintf('%s/hsc/stackdathsc_%s_masklim_bk',...
                loaddir,name),'stackdathscbk');  
        else
            load(sprintf('%s/hsc/stackdathsc_%s',...
                loaddir,name),'stackdathsc');
            load(sprintf('%s/hsc/stackdathsc_%s_bk',...
                loaddir,name),'stackdathscbk');
        end
        
        stackdat = stackdathsc(im).stackdat;
        bkallcb_arr = ones([numel(stackdathscbk),numel(r_arr)]).*nan;
        bkallps_arr = ones([numel(stackdathscbk),numel(r_arr)]).*nan;
        for iter=1:numel(stackdathscbk)
            stackdatbk = stackdathscbk(iter).stackdat(im);
            sp  = find(stackdatbk.hitmapg_arr~=0);
            bkallcb_arr(iter,sp) = stackdatbk.profcbg_arr(sp);
            bkallps_arr(iter,sp) = stackdatbk.profpsg_arr(sp);
        end
        bkavgcb_arr = nanmean(bkallcb_arr);
        sp = find(bkavgcb_arr==bkavgcb_arr);
        bkavgcb_arr = spline(r_arr(sp),bkavgcb_arr(sp),r_arr);
        bkavgcb_err = nanstd(bkallcb_arr);
        sp = find(bkavgcb_err==bkavgcb_err);
        bkavgcb_err = spline(r_arr(sp),bkavgcb_err(sp),r_arr);
        bkavgps_arr = nanmean(bkallps_arr);
        sp = find(bkavgps_arr==bkavgps_arr);
        bkavgps_arr = spline(r_arr(sp),bkavgps_arr(sp),r_arr);
        bkavgps_err = nanstd(bkallps_arr);
        sp = find(bkavgps_err==bkavgps_err);
        bkavgps_err = spline(r_arr(sp),bkavgps_err(sp),r_arr);

        profcbg = stackdat.all.profcbg - bkavgcb_arr;
        profcbg_err = sqrt(stackdat.errjack.profcbg.^2 + bkavgcb_err.^2);
        profpsg = stackdat.all.profpsg - bkavgps_arr;
        profpsg_err = sqrt(stackdat.errjack.profpsg.^2 + bkavgps_err.^2);
        excb_all(hsc_idx+1,:) = profcbg - psf_arr.*profcbg(1);
        exps_all(hsc_idx+1,:) = profpsg - psf_arr.*profpsg(1);
        excb_err_all(hsc_idx+1,:) = profcbg_err;
        exps_err_all(hsc_idx+1,:) = profpsg_err;
    end
    excb_all = sum(excb_all./(excb_err_all.^2))./sum(1./(excb_err_all.^2));
    exps_all = sum(exps_all./(exps_err_all.^2))./sum(1./(excb_err_all.^2));
    excb_err_all = sqrt(1./sum(1./(excb_err_all.^2)));
    exps_err_all = sqrt(1./sum(1./(exps_err_all.^2)));
    excess(2).cb(im).ihl = excb_all;
    excess(2).cb(im).err = excb_err_all;
    excess(2).ps(im).ihl = exps_all;
    excess(2).ps(im).err = exps_err_all;
end

%%
%%% plot the stacking profile %%%

for ifield = 4:8
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
setwinsize(gcf,1200,750)
for im=1:4
    stackdat = stackdatall(im).stackdat;
    
    m_min = stackdat.m_min;
    m_max = stackdat.m_max;
    counts = stackdat.all.counts;
    countg = stackdat.all.countg;
    r_arr = stackdat.r_arr;

    % profscb = stackdat.all.profcbs - cbbk;
    % profscb_err = stackdat.errjack.profcbs;
    % profsps = stackdat.all.profpss - psbk;
    % profsps_err = stackdat.errjack.profpss;
    profscb = stackdat.norm.profcbpsf;
    profscb_err = stackdat.norm.profcbpsf_err;
    profsps = stackdat.norm.profpspsf;
    profsps_err = stackdat.norm.profcbpsf_err;
    
    cbmean = mean(stackdat.bk.profcbg(r_arr > 100));
    psmean = mean(stackdat.bk.profpsg(r_arr > 100));
    profgcb = stackdat.all.profcbg - cbmean;
    profgcb_err = stackdat.errjack.profcbg;
    profgps = stackdat.all.profpsg - psmean;
    profgps_err = stackdat.errjack.profpsg;

    profgcbbk = stackdat.bk.profcbg - cbmean;
    profgcb_errbk = stackdat.bk.profcbg_err;
    profgpsbk = stackdat.bk.profpsg - psmean;
    profgps_errbk = stackdat.bk.profpsg_err;

    subplot(3,4,im)
    semilogx(r_arr.*0.98,profscb,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgcb,'b','markersize',10);
    plot(r_arr.*0.98,profsps,'m','markersize',10);
    plot(r_arr.*1.02,profgps,'c','markersize',10);
    h=legend({'CIBER PSF','CIBER galaxies',...
        'PanSTARRS stars','PanSTARRS galaxies'},'Location','southwest');
    set(h,'fontsize',7)
    legend boxoff
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    errorbar(r_arr.*0.98,profsps,profsps_err,'m.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'c.','markersize',10);
    xlim([4e-1,50])
    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
    title({strcat(num2str(m_min),'<m<',num2str(m_max)),...
        strcat(num2str(counts),' stars, ',num2str(countg), ' galaxies')},...
        'fontsize',15)
            
    subplot(3,4,im + 4)
    semilogx(r_arr.*0.98,profscb,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgcb,'b','markersize',10);
    plot(r_arr,profgcbbk,'k','markersize',10);
    h=legend({'CIBER PSF','CIBER galaxies','background (gals)'},...
        'Location','northeast');
    set(h,'fontsize',7)
    legend boxoff
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    errorbar(r_arr,profgcbbk,profgcb_errbk,'k.','markersize',10);
    xlim([4e-1,1.1e3])
    ylim([-5,10])
    if masklim
        ylim([profgcbbk(end)-40*profgcb_errbk(end),...
            profgcbbk(end)+100*profgcb_errbk(end)])
    end
    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)

    subplot(3,4,im + 8)
    semilogx(r_arr.*0.98,profsps,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgps,'b','markersize',10);
    plot(r_arr,profgpsbk,'k','markersize',10);
    h=legend({'PanSTARRS stars','PanSTARRS galaxies','background (gals)'},...
        'Location','northeast');
    set(h,'fontsize',7)
    legend boxoff
    errorbar(r_arr.*0.98,profsps,profsps_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'b.','markersize',10);
    errorbar(r_arr,profgpsbk,profgps_errbk,'k.','markersize',10);
    xlim([4e-1,1.1e3])
    ylim([-0.01,0.04])
    if masklim
        ylim([profgpsbk(end)-40*profgps_errbk(end),...
            profgpsbk(end)+100*profgps_errbk(end)])
    end
    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end

if masklim
    suptitle(strcat(dt.name,' (mask to mag bin max)'));
else
    suptitle(strcat(dt.name,' (mask all PanSTARRS sources)'));
end

if masklim
    savename = sprintf('%s/%s_stackprof_masklim',...
        pltsavedir,dt.name);
else
    savename = sprintf('%s/%s_stackprof',...
        pltsavedir,dt.name);    
end
if savefig
    print(savename,'-dpng');close
end

end

%%
%%% plot the field-averaged stacking profile %%%

figure
setwinsize(gcf,1200,750)
for im=1:4
    profscb = 0;
    profscb_err = 0;
    profgcb = 0;
    profgcb_err = 0;
    profsps = 0;
    profsps_err = 0;
    profgps = 0;
    profgps_err = 0;
    for ifield = 4:8
        dt=get_dark_times(flight,inst,ifield);
        loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
        if masklim
            load(sprintf('%s/stackdat_%s_masklim',...
                    loaddir,dt.name),'stackdatall');        
        else
            load(sprintf('%s/stackdat_%s',...
                    loaddir,dt.name),'stackdatall');
        end
        stackdat = stackdatall(im).stackdat;
        
        profscb = profscb + ...
            stackdat.norm.profcbpsf./(stackdat.norm.profcbpsf_err.^2);
        profscb_err = profscb_err + ...
            1./(stackdat.norm.profcbpsf_err.^2);
        profgcb = profgcb + ...
            stackdat.norm.profcbg./(stackdat.norm.profcbg_err.^2);
        profgcb_err = profgcb_err + ...
            1./(stackdat.norm.profcbg_err.^2);
        profsps = profsps + ...
            stackdat.norm.profpspsf./(stackdat.norm.profpspsf_err.^2);
        profsps_err = profsps_err + ...
            1./(stackdat.norm.profpspsf_err.^2);
        profgps = profgps + ...
            stackdat.norm.profpsg./(stackdat.norm.profpsg_err.^2);
        profgps_err = profgps_err + ...
            1./(stackdat.norm.profpsg_err.^2);                    
    end
    profscb = profscb./profscb_err;
    profscb_err = sqrt(1./profscb_err);
    profgcb = profgcb./profgcb_err;
    profgcb_err = sqrt(1./profgcb_err);
    profsps = profsps./profsps_err;
    profsps_err = sqrt(1./profsps_err);
    profgps = profgps./profgps_err;
    profgps_err = sqrt(1./profgps_err);
    
    m_min = stackdat.m_min;
    m_max = stackdat.m_max;
    r_arr = stackdat.r_arr;

    subplot(3,4,im)
    semilogx(r_arr.*0.98,profscb,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgcb,'b','markersize',10);
    plot(r_arr.*0.98,profsps,'m','markersize',10);
    plot(r_arr.*1.02,profgps,'c','markersize',10);
    h=legend({'CIBER stars','CIBER galaxies',...
        'PanSTARRS stars','PanSTARRS galaxies'},'Location','southwest');
    set(h,'fontsize',7)
    legend boxoff
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    errorbar(r_arr.*0.98,profsps,profsps_err,'m.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'c.','markersize',10);
    xlim([4e-1,50])
    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
    title({strcat(num2str(m_min),'<m<',num2str(m_max))},...
        'fontsize',15)
            
    subplot(3,4,im + 4)
    semilogx(r_arr.*0.98,profscb,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgcb,'b','markersize',10);
    h=legend({'CIBER stars','CIBER galaxies'},...
        'Location','northeast');
    set(h,'fontsize',7)
    legend boxoff
    errorbar(r_arr.*0.98,profscb,profscb_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgcb,profgcb_err,'b.','markersize',10);
    xlim([4e-1,1.1e3])
    ylim([-2,10])
    if masklim
        ylim([profscb(end)-40*profscb_err(end),...
            profscb(end)+100*profscb_err(end)])
    end

    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)

    subplot(3,4,im + 8)
    semilogx(r_arr.*0.98,profsps,'r','markersize',10);hold on
    plot(r_arr.*1.02,profgps,'b','markersize',10);
    h=legend({'PanSTARRS stars','PanSTARRS galaxies'},...
        'Location','northeast');
    set(h,'fontsize',7)
    legend boxoff
    errorbar(r_arr.*0.98,profsps,profsps_err,'r.','markersize',10);
    errorbar(r_arr.*1.02,profgps,profgps_err,'b.','markersize',10);
    xlim([4e-1,1.1e3])
    ylim([-0.01,0.04])
    if masklim
        ylim([profsps(end)-400*profsps_err(end),...
            profsps(end)+1000*profsps_err(end)])
    end
    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end

if masklim
    suptitle(strcat('field-averaged profile (mask to mag bin max)'));
else
    suptitle(strcat('field-averaged profile (mask all PanSTARRS sources)'));
end

if masklim
    savename = sprintf('%s/stackprof_masklim',...
        pltsavedir);
else
    savename = sprintf('%s/stackprof',...
        pltsavedir);    
end

if savefig
    print(savename,'-dpng');close
end
%%
%%% plot the excess profile in each field %%%

for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

figure
setwinsize(gcf,1400,700)
for im=1:4
    m_min = excess(ifield).mag(im).m_min;
    m_max = excess(ifield).mag(im).m_max;
    counts = excess(ifield).count(im).counts;
    countg = excess(ifield).count(im).countg;
    r_arr = excess(ifield).r_arr;

    ecb = excess(ifield).cb(im).ihl;
    ecb_err = excess(ifield).cb(im).err;
    eps = excess(ifield).ps(im).ihl;
    eps_err = excess(ifield).ps(im).err;
    ehsc = excess(2).cb(im).ihl;
    ehsc_err = excess(2).cb(im).err;
    
    subplot(2,4,im)
    loglog(r_arr,ecb,'k.','markersize',10);hold on
    loglog(r_arr.*1.02,eps,'b.','markersize',10);hold on
    loglog(r_arr.*0.98,ehsc,'r.','markersize',10);hold on
    if im==4
        h=legend({'CIBER Excess','PanSTARRS Excess', 'HSC Excess',},...
            'Location','northeast');
        set(h,'fontsize',10)
    end
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
    errorbar(r_arr,ecb,ecb_err_low,ecb_err,'k.','markersize',10);
    loglog(r_arr,-ecb,'ko','markersize',5);
    ecb_err_low = ecb_err;
    for ii=1:numel(ecb)
        if ecb(ii)<0 & ecb_err(ii) > abs(ecb(ii))
            ecb_err_low(ii) = -ecb(ii) - 1e-8;
        end
    end
    errorbar(r_arr,-ecb,ecb_err_low,ecb_err,'ko','markersize',5);
    
    ehsc_err_low = ehsc_err;
    for ii=1:numel(eps)
        if ehsc(ii)>0 & ehsc_err(ii) > ehsc(ii)
            ehsc_err_low(ii) = ehsc(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*0.98,ehsc,ehsc_err_low,ehsc_err,'r.','markersize',10);
    loglog(r_arr.*0.98,-ehsc,'ro','markersize',5);
    ehsc_err_low = ehsc_err;
    for ii=1:numel(ehsc)
        if ehsc(ii)<0 & ehsc_err(ii) > abs(ehsc(ii))
            ehsc_err_low(ii) = -ehsc(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*0.98,-ehsc,ehsc_err_low,ehsc_err,'ro','markersize',5);
    
    
    xlim([4e-1,1.1e3])
    ylim([1e-4,1e4])
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
      ' (',num2str(counts),' stars, ',num2str(countg), ' galaxies)'),'fontsize',10)
    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
    
    subplot(2,4,im+4)
    errorbar(r_arr, ecb, excess(ifield).cb(im).errdf,'r.','markersize',10);hold on
    errorbar(r_arr, ecb, ecb_err,'k.','markersize',10);
    sp = find(r_arr > 100);
    v = excess(ifield).cb(im).ihl100;
    e = excess(ifield).cb(im).err100;
    plot(r_arr(sp),v*ones(size(sp)),'color',[0.2,0.4,0.2],'linewidth',2);
    if im==4
        h=legend({'CIBER Excess (w/ fluc.)','CIBER Excess',...
            'CIBER Excess > 100 arcsec',},'Location','northeast');
        set(h,'fontsize',10)
    end
    fill([r_arr(sp),flip(r_arr(sp))],...
    [(v+e) * ones(size(sp)),(v-e) * ones(size(sp))],...
        [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
    set(gca, 'XScale', 'log');
    hline(0,'k:');
    xlim([4e-1,1.1e3])
    ylim([-2,8])
    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)

end

if masklim
    suptitle(strcat(dt.name,' (mask to mag bin max)'));
else
    suptitle(strcat(dt.name,' (mask all PanSTARRS sources)'));
end

if masklim
    savename = sprintf('%s/%s_excessprof_masklim',...
        pltsavedir,dt.name);
else
    savename = sprintf('%s/%s_excessprof',...
        pltsavedir,dt.name);    
end

if savefig
    print(savename,'-dpng');close
end
end

%%
%%% plot the excess profile >100 arcsec in each field %%%

figure
setwinsize(gcf,1400,300)
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
for im=1:4
    subplot(1,4,im)
    name = {};
    for ifield = 4:8
        dt=get_dark_times(flight,inst,ifield);
        name{end+1} = dt.name;
        m_min = excess(ifield).mag(im).m_min;
        m_max = excess(ifield).mag(im).m_max;
        r_arr = excess(ifield).r_arr;

        sp = find(r_arr > 100);
        v = excess(ifield).cb(im).ihl100;
        e = excess(ifield).cb(im).err100;
        errorbar(ifield-3,v,e,'k.','markersize',10);hold on
    end
    plot([0,6],[0,0],'b--')
    ylabel('mean excess > 100 arcsec [nW/m^2/sr]', 'fontsize',10);
    xlim([0.5,5.5])
    set(gca,'xTick',1:5);
    set(gca, 'xTickLabels', name); 
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);

end


if masklim
    savename = sprintf('%s/excess100_masklim',...
        pltsavedir);
else
    savename = sprintf('%s/excess100',...
        pltsavedir);    
end
if savefig
    print(savename,'-dpng');close
end
%% plot all field excess in a plot

figure
setwinsize(gcf,1500,600)

for im=1:4
ihldat = zeros([5,numel(r_arr)]);
errdat = zeros([5,numel(r_arr)]);

for ifield = 4:8
    dt=get_dark_times(flight,inst,ifield);

    m_min = excess(ifield).mag(im).m_min;
    m_max = excess(ifield).mag(im).m_max;
    r_arr = excess(ifield).r_arr;

    ihl = excess(ifield).cb(im).ihl;
    ihl_err = excess(ifield).cb(im).err;
    ihldat(ifield-3,:) = ihl;
    errdat(ifield-3,:) = ihl_err;
    
    off = 1 + (ifield - 6).*0.02;
    subplot(2,4,im)
    loglog(r_arr.*off,ihl,'.','color',...
        get_color(ifield-3),'markersize',10,'DisplayName',dt.name);hold on
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);

    
end
if im==4
    h=legend('show','Location','northeast');
    set(h,'fontsize',10)
    legend boxoff
end

ihlmid = median(ihldat);
errmid = median(errdat);
for ifield = 4:8
    ihl = excess(ifield).cb(im).ihl;
    ihl_err = excess(ifield).cb(im).err; 
    off = 1 + (ifield - 6).*0.02;
    subplot(2,4,im)
    
    ihl_err_low = ihl_err;
    for ii=1:numel(ihl)
        if ihl(ii)>0 & ihl_err(ii) > ihl(ii)
            ihl_err_low(ii) = ihl(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*off,ihl,ihl_err_low,ihl_err,'.','color',...
        get_color(ifield-3),'markersize',10);hold on
    
    loglog(r_arr.*off,-ihl,'o','color',...
        get_color(ifield-3),'markersize',5);hold on
    ihl_err_low = ihl_err;
    for ii=1:numel(ihl)
        if ihl(ii)<0 & ihl_err(ii) > abs(ihl(ii))
            ihl_err_low(ii) = -ihl(ii) - 1e-8;
        end
    end    
    errorbar(r_arr.*off,-ihl,ihl_err_low,ihl_err,'o','color',...
        get_color(ifield-3),'markersize',5);hold on

    xlim([4e-1,1.1e3])
    ylim([1e-4,1e3])
    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end

for ifield = 4:8

    ihl = excess(ifield).cb(im).ihl;
    ihl_err = excess(ifield).cb(im).err;

    off = 1 + (ifield - 6).*0.02;
    subplot(2,4,im + 4)
    semilogx(r_arr.*off,(ihl - ihlmid)./errmid,'.','color',...
        get_color(ifield-3),'markersize',10);hold on
    errorbar(r_arr.*off,(ihl - ihlmid)./errmid,ihl_err./errmid,'.','color',...
        get_color(ifield-3),'markersize',10);
    xlim([4e-1,1.1e3])
    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('scaled excess', 'fontsize',15)

end
end

if masklim
    savename = sprintf('%s/excess_all_masklim',...
        pltsavedir);
else
    savename = sprintf('%s/excess_all',...
        pltsavedir);    
end

if savefig
    print(savename,'-dpng');close
end
%% plot field-average excess

figure
setwinsize(gcf,1400,700)

for im=1:4
    m_min = excess(ifield).mag(im).m_min;
    m_max = excess(ifield).mag(im).m_max;
    r_arr = excess(ifield).r_arr;
    ecb = excess(1).cb(im).ihl;
    ecb_err = excess(1).cb(im).err;
    eps = excess(1).ps(im).ihl;
    eps_err = excess(1).ps(im).err;
    ehsc = excess(2).cb(im).ihl;
    ehsc_err = excess(2).cb(im).err;

    subplot(2,4,im)
    loglog(r_arr,ecb,'k.','markersize',10);hold on
    loglog(r_arr.*1.05,eps,'b.','markersize',10);hold on
    loglog(r_arr.*0.95,ehsc,'r.','markersize',10);hold on
    if im==4
        h=legend({'CIBER Excess','PanSTARRS Excess', 'HSC Excess',},...
            'Location','northeast');
        set(h,'fontsize',10)
    end
    
    eps_err_low = eps_err;
    for ii=1:numel(eps)
        if eps(ii)>0 & eps_err(ii) > eps(ii)
            eps_err_low(ii) = eps(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*1.05,eps,eps_err_low,eps_err,'b.','markersize',10);
    loglog(r_arr.*1.05,-eps,'bo','markersize',5);hold on
    eps_err_low = eps_err;
    for ii=1:numel(eps)
        if eps(ii)<0 & eps_err(ii) > abs(eps(ii))
            eps_err_low(ii) = -eps(ii) - 1e-8;
        end
    end    
    errorbar(r_arr.*1.05,-eps,eps_err_low,eps_err,'bo','markersize',5);
    
    ecb_err_low = ecb_err;
    for ii=1:numel(eps)
        if ecb(ii)>0 & ecb_err(ii) > ecb(ii)
            ecb_err_low(ii) = ecb(ii) - 1e-8;
        end
    end
    errorbar(r_arr,ecb,ecb_err_low,ecb_err,'k.','markersize',10);
    loglog(r_arr,-ecb,'ko','markersize',5);
    ecb_err_low = ecb_err;
    for ii=1:numel(ecb)
        if ecb(ii)<0 & ecb_err(ii) > abs(ecb(ii))
            ecb_err_low(ii) = -ecb(ii) - 1e-8;
        end
    end
    errorbar(r_arr,-ecb,ecb_err_low,ecb_err,'ko','markersize',5);
    
    ehsc_err_low = ehsc_err;
    for ii=1:numel(eps)
        if ehsc(ii)>0 & ehsc_err(ii) > ehsc(ii)
            ehsc_err_low(ii) = ehsc(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*0.95,ehsc,ehsc_err_low,ehsc_err,'r.','markersize',10);
    loglog(r_arr.*0.95,-ehsc,'ro','markersize',5);
    ehsc_err_low = ehsc_err;
    for ii=1:numel(ehsc)
        if ehsc(ii)<0 & ehsc_err(ii) > abs(ehsc(ii))
            ehsc_err_low(ii) = -ehsc(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*0.95,-ehsc,ehsc_err_low,ehsc_err,'ro','markersize',5);
    ylim([1e-4,1e3])
    xlim([4e-1,1.1e3])
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);
    xlabel('r [arcsec]', 'fontsize',15);
    ylabel('I [nW/m^2/sr]', 'fontsize',15);

    subplot(2,4,im+4)
    errorbar(r_arr, ecb, ecb_err,'k.','markersize',10);hold on
    sp = find(r_arr > 100);
    v = excess(1).cb(im).ihl100;
    e = excess(1).cb(im).err100;
    plot(r_arr(sp),v*ones(size(sp)),'color',[0.2,0.4,0.2],'linewidth',2);
    if im==4
        h=legend({'CIBER Excess','CIBER Excess > 100 arcsec',},...
            'Location','northeast');
        set(h,'fontsize',10)
    end
    fill([r_arr(sp),flip(r_arr(sp))],...
    [(v+e) * ones(size(sp)),(v-e) * ones(size(sp))],...
        [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');
    set(gca, 'XScale', 'log');
    hline(0,'k:');
    xlim([4e-1,1.1e3])
    ylim([-1,3])
    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)
end
if masklim
    savename = sprintf('%s/excess_avg_masklim',...
        pltsavedir);
else
    savename = sprintf('%s/excess_avg',...
        pltsavedir);    
end

if savefig
    print(savename,'-dpng');close
end
%% plot field-average integrated excess
figure
setwinsize(gcf,1400,700)
r_arr = excess(4).r_arr;
dx=1200;
nbins=25;
profile = radial_prof(ones(2*dx+1),ones(2*dx+1),dx+1,dx+1,1,nbins);
rbinedges_arr= profile.binedges.*0.7;
rbinedges_arr = rbinedges_arr.*pi./180/3600;% [rad]
dA = pi*rbinedges_arr(2:end).^2-pi*rbinedges_arr(1:end-1).^2; % [sr]
for im=1:4
    m_min = excess(ifield).mag(im).m_min;
    m_max = excess(ifield).mag(im).m_max;
    
    ecb = cumsum(excess(1).cb(im).ihl.*dA).*1e9;% [W/m^2]
    ecb_err = sqrt(cumsum(excess(1).cb(im).err.^2.*dA.^2)).*1e9;
    eps = cumsum(excess(1).ps(im).ihl.*dA).*1e9;
    eps_err = sqrt(cumsum(excess(1).ps(im).err.^2.*dA.^2)).*1e9;
    ehsc = cumsum(excess(2).cb(im).ihl.*dA).*1e9;
    ehsc_err = sqrt(cumsum(excess(2).cb(im).err.^2.*dA.^2)).*1e9;
    src = cumsum(excess(1).cb(im).src.*dA).*1e9;
    src_err = sqrt(cumsum(excess(1).cb(im).src_err.^2.*dA.^2)).*1e9;
    
    subplot(2,4,im)
    loglog(r_arr,ecb,'k.','markersize',10);hold on
    loglog(r_arr.*1.05,eps,'b.','markersize',10);
    loglog(r_arr.*0.95,ehsc,'r.','markersize',10);
    loglog(r_arr,src,'m.','markersize',10);
    if im==4
        h=legend({'CIBER Excess','PanSTARRS Excess','HSC Excess','Source'},...
            'Location','northwest');
        set(h,'fontsize',8)
    end

    ecb_err_low = ecb_err;
    for ii=1:numel(ecb)
        if ecb(ii)>0 & ecb_err(ii) > ecb(ii)
            ecb_err_low(ii) = ecb(ii) - 1e-8;
        end
    end
    errorbar(r_arr,ecb,ecb_err_low,ecb_err,'k.','markersize',10);
    loglog(r_arr,-ecb,'ko','markersize',5);hold on
    ecb_err_low = ecb_err;
    for ii=1:numel(ecb)
        if ecb(ii)<0 & ecb_err(ii) > abs(ecb(ii))
            ecb_err_low(ii) = -ecb(ii) - 1e-8;
        end
    end    
    errorbar(r_arr,-ecb,ecb_err_low,ecb_err,'ko','markersize',5);
    
    eps_err_low = eps_err;
    for ii=1:numel(eps)
        if eps(ii)>0 & eps_err(ii) > eps(ii)
            eps_err_low(ii) = eps(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*1.05,eps,eps_err_low,eps_err,'b.','markersize',10);
    loglog(r_arr.*1.05,-eps,'bo','markersize',5);hold on
    eps_err_low = eps_err;
    for ii=1:numel(eps)
        if eps(ii)<0 & eps_err(ii) > abs(eps(ii))
            eps_err_low(ii) = -eps(ii) - 1e-8;
        end
    end    
    errorbar(r_arr.*1.05,-eps,eps_err_low,eps_err,'bo','markersize',5);

    ehsc_err_low = ehsc_err;
    for ii=1:numel(eps)
        if ehsc(ii)>0 & ehsc_err(ii) > ehsc(ii)
            ehsc_err_low(ii) = ehsc(ii) - 1e-8;
        end
    end
    errorbar(r_arr.*0.95,ehsc,ehsc_err_low,ehsc_err,'r.','markersize',10);
    loglog(r_arr.*0.95,-ehsc,'ro','markersize',5);hold on
    ehsc_err_low = ehsc_err;
    for ii=1:numel(ehsc)
        if ehsc(ii)<0 & ehsc_err(ii) > abs(ehsc(ii))
            ehsc_err_low(ii) = -ehsc(ii) - 1e-8;
        end
    end    
    errorbar(r_arr.*0.95,-ehsc,ehsc_err_low,ehsc_err,'ro','markersize',5);

    
    src_err_low = src_err;
    for ii=1:numel(src)
        if src(ii)>0 & src_err(ii) > src(ii)
            src_err_low(ii) = src(ii) - 1e-8;
        end
    end
    errorbar(r_arr,src,src_err_low,src_err,'m.','markersize',10);
    loglog(r_arr,-src,'mo','markersize',5);hold on
    src_err_low = src_err;
    for ii=1:numel(src)
        if src(ii)<0 & src_err(ii) > abs(src(ii))
            src_err_low(ii) = -src(ii) - 1e-8;
        end
    end    
    errorbar(r_arr,-src,src_err_low,src_err,'mo','markersize',5);
    
    ylim([1e-2,1e5])
    xlim([4e-1,1.1e3])
    xlabel('r [arcsec]', 'fontsize',15)
    ylabel('$\nu F_\nu\,{\rm[W/m^2]}$',...
    'interpreter','latex', 'fontsize',15);
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);

    subplot(2,4,im+4)
    loglog(1e-8, 1e-8,'k.','markersize',10);hold on
    plot(1e-8, 1e-8,'b.','markersize',10);
    plot(1e-8, 1e-8,'r.','markersize',10);
    if im==4
    h=legend({'CIBER Excess','PanSTARRS Excess','HSC Excess'},...
        'Location','northwest');
    set(h,'fontsize',8)
    end
    
    sp = 1:18;
    
    r = r_arr(sp);    
    v = ecb(sp)./src(sp);
    e = abs(ecb./src).*sqrt((ecb_err./ecb).^2+(src_err./src).^2);
    e = e(sp);
    loglog(r, v,'k.','markersize',10);hold on
    errorbar(r, v, e,'k.','markersize',10);
    elow = e;
    elow(v>=0 & v<e) = v(v>=0 & v<e)-1e-8;
    errorbar(r, v, elow, e,'k.','markersize',10);
    errorbar(r, -v, e,'ko','markersize',5);
    elow(v<0 & -v<e) = -v(v<0 & -v<e)-1e-8;
    errorbar(r, -v, elow, e,'ko','markersize',5);
    
    r = r_arr(sp).*1.05;
    v = eps(sp)./src(sp);
    e = abs(eps./src).*sqrt((eps_err./eps).^2+(src_err./src).^2);
    e = e(sp);
    loglog(r, v,'b.','markersize',10);hold on
    errorbar(r, v, e,'b.','markersize',10);
    elow = e;
    elow(v>=0 & v<e) = v(v>=0 & v<e)-1e-8;
    errorbar(r, v, elow, e,'b.','markersize',10);
    errorbar(r, -v, e,'bo','markersize',5);
    elow(v<0 & -v<e) = -v(v<0 & -v<e)-1e-8;
    errorbar(r, -v, elow, e,'bo','markersize',5);
    
    r = r_arr(sp).*0.95;
    v = ehsc(sp)./src(sp);
    e = abs(ehsc./src).*sqrt((ehsc_err./ehsc).^2+(src_err./src).^2);
    e = e(sp);
    loglog(r, v,'r.','markersize',10);hold on
    errorbar(r, v, e,'r.','markersize',10);
    elow = e;
    elow(v>=0 & v<e) = v(v>=0 & v<e)-1e-8;
    errorbar(r, v, elow, e,'r.','markersize',10);
    errorbar(r, -v, e,'ro','markersize',5);
    elow(v<0 & -v<e) = -v(v<0 & -v<e)-1e-8;
    errorbar(r, -v, elow, e,'ro','markersize',5);

    hline(0.1,'k--');
    hline(0.3,'k--');
    hline(0.5,'k--');
    ylim([1e-3,1e1])
    xlim([4e-1,1.1e3])
    xlabel('r [arcsec]', 'fontsize',15);
    ylabel('$\nu F_\nu$\,(excess)\,/\,$\nu F_\nu$\,(source)',...
    'interpreter','latex', 'fontsize',15);

end

if masklim
    savename = sprintf('%s/excessint_masklim',...
        pltsavedir);
else
    savename = sprintf('%s/excessint',...
        pltsavedir);    
end
if savefig
    print(savename,'-dpng');close
end
%% total EBL intensity from excess from each mag bin
figure
r_arr = excess(4).r_arr;
sp10 = find(abs(r_arr-10)==min(abs(r_arr-10)));
sp100 = find(abs(r_arr-100)==min(abs(r_arr-100)));
dx=1200;
nbins=25;
profile = radial_prof(ones(2*dx+1),ones(2*dx+1),dx+1,dx+1,1,nbins);
rbinedges_arr= profile.binedges.*0.7;
rbinedges_arr = rbinedges_arr.*pi./180/3600;% [rad]
dA = pi*rbinedges_arr(2:end).^2-pi*rbinedges_arr(1:end-1).^2; % [sr]
I10_cb = zeros([1,4]);
I100_cb = zeros([1,4]);
e10_cb = zeros([1,4]);
e100_cb = zeros([1,4]);
I10_hsc = zeros([1,4]);
I100_hsc = zeros([1,4]);
e10_hsc = zeros([1,4]);
e100_hsc = zeros([1,4]);

for im=1:4
    m_min = excess(ifield).mag(im).m_min;
    m_max = excess(ifield).mag(im).m_max;    
    if inst==1
        nu = 3e8/1.1e-6; % [Hz]
    else
        nu = 3e8/1.6e-6; % [Hz]
    end
    nuFnu = nu*3631*10^(-(m_min+m_max)/2/2.5)*1e-17; % [nW/m^2]
    nuInu = nuFnu*IGLcounts_helgason(inst,(m_min+m_max)/2)*(180/pi)^2;% [nW/m^2/sr]

    ecb = cumsum(excess(1).cb(im).ihl.*dA).*1e9;% [W/m^2]
    ecb_err = sqrt(cumsum(excess(1).cb(im).err.^2.*dA.^2)).*1e9;
    ehsc = cumsum(excess(2).cb(im).ihl.*dA).*1e9;
    ehsc_err = sqrt(cumsum(excess(2).cb(im).err.^2.*dA.^2)).*1e9;
    src = cumsum(excess(1).cb(im).src.*dA).*1e9;
    src_err = sqrt(cumsum(excess(1).cb(im).src_err.^2.*dA.^2)).*1e9;
    
    I10_cb(im) = ecb(sp10)./src(sp10).*nuInu;
    I100_cb(im) = ecb(sp100)./src(sp100).*nuInu;
    e=abs(ecb./src).*sqrt((ecb_err./ecb).^2+(src_err./src).^2);
    e10_cb(im) = e(sp10).*nuInu;
    e100_cb(im) = e(sp100).*nuInu;

    I10_hsc(im) = ehsc(sp10)./src(sp10).*nuInu;
    I100_hsc(im) = ehsc(sp100)./src(sp100).*nuInu;
    e=abs(ehsc./src).*sqrt((ehsc_err./ehsc).^2+(src_err./src).^2);
    e10_hsc(im) = e(sp10).*nuInu;
    e100_hsc(im) = e(sp100).*nuInu;
    
end

errorbar(16.52:19.52,I10_cb,e10_cb,'ko-','markersize',5);hold on
errorbar(16.49:19.49,I10_hsc,e10_hsc,'ro-','markersize',5);hold on
errorbar(16.5:19.5,I100_cb,e100_cb,'k.-','markersize',10);
errorbar(16.51:19.51,I100_hsc,e100_hsc,'r.-','markersize',10);
h=legend({'CIBER<10 arcsec','HSC<10 arcsec','CIBER<100 arcsec',...
    'HSC<100 arcsec'},'Location','northwest');
set(h,'fontsize',12)

ylim([0,3])
xlim([16,20])
xlabel('mag', 'fontsize',15)
ylabel('I[nW/m^2/sr]', 'fontsize',15)

if masklim
    savename = sprintf('%s/IEBL_masklim',...
        pltsavedir);
else
    savename = sprintf('%s/IEBL',...
        pltsavedir);    
end
if savefig
    print(savename,'-dpng');close
end

%%
return