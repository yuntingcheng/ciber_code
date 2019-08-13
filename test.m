flight=40030;
mypaths=get_paths(flight);
inst=1;
masklim=true;
m_min_arr = 16:19;
m_max_arr = m_min_arr + 1;

%%
% [psf_arr,~,~] =  PSF_stacked_profile(flight,inst,8);
%% process data
profcb_all = zeros([4,25]);
profcb_all_err = zeros([4,25]);
profps_all = zeros([4,25]);
profps_all_err = zeros([4,25]);
ecb_all = zeros([4,25]);
eps_all = zeros([4,25]);
for hsc_idx=0:11

    [name,~] = HSC_fields_info(hsc_idx);
    loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

    if masklim
        load(sprintf('%s/stackdathsc_%s_masklim',...
            loaddir,name),'stackdathsc');  
        load(sprintf('%s/stackdathsc_%s_masklim_bk',...
            loaddir,name),'stackdathscbk');  

    else
        load(sprintf('%s/stackdathsc_%s',...
            loaddir,name),'stackdathsc');
        load(sprintf('%s/stackdathsc_%s_bk',...
            loaddir,name),'stackdathscbk');
    end

    for im=1:4 
        stackdat = stackdathsc(im).stackdat;
        r_arr = stackdat.r_arr;
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
        
        profcb_all(im,:) = profcb_all(im,:) + profcbg./profcbg_err.^2;
        profcb_all_err(im,:) = profcb_all_err(im,:) + 1./profcbg_err.^2;
        profps_all(im,:) = profps_all(im,:) + profpsg./profpsg_err.^2;
        profps_all_err(im,:) = profps_all_err(im,:) + 1./profpsg_err.^2;
        
        dat(im).r_arr = r_arr;
        dat(im).psf_arr = psf_arr;
        dat(im).Ng = stackdat.all.countg;
        dat(im).profcbg = profcbg;
        dat(im).errcb = profcbg_err;
        dat(im).profpsg = profpsg;
        dat(im).errps = profpsg_err;
        dat(im).excesscb = profcbg - psf_arr.*profcbg(1);
        dat(im).excessps = profpsg - psf_arr.*profpsg(1);
        
        ecb_all(im,:) = ecb_all(im,:) +dat(im).excesscb./profcbg_err.^2;
        eps_all(im,:) = eps_all(im,:) +dat(im).excessps./profpsg_err.^2;
    end
    datall(hsc_idx+1).dat = dat;
end

profcb_all = profcb_all./profcb_all_err;
ecb_all = ecb_all./profcb_all_err;
profcb_all_err = sqrt(1./profcb_all_err);
profps_all = profps_all./profps_all_err;
eps_all = eps_all./profps_all_err;
profps_all_err = sqrt(1./profps_all_err);
%%

dat = datall(3).dat;
[name,~] = HSC_fields_info(3);
figure
setwinsize(gcf,1200,750)

for im=1:4
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    r_arr = dat(im).r_arr;
    psf_arr = dat(im).psf_arr;
    Ng = dat(im).Ng;
    
    subplot(3,4,im)
    semilogx(r_arr,psf_arr.*dat(im).profcbg(1),'k');hold on
    errorbar(r_arr,dat(im).profcbg,dat(im).errcb,'r.');
    errorbar(r_arr,dat(im).profpsg,dat(im).errps,'b.');
    xlim([4e-1,50]);
    if im==1
        ylabel('I [nW/m^2/sr]', 'fontsize',15);
    end
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
        '(',num2str(Ng), ' galaxies)'),'fontsize',15);

    if im==4
        h=legend({'PSF','HSC all','HSC m<20'},'Location','northeast');
        set(h,'fontsize',7)
    end
    subplot(3,4,im+4)
    semilogx(r_arr,psf_arr.*dat(im).profcbg(1),'k');hold on
    errorbar(r_arr.*0.99,dat(im).profcbg,dat(im).errcb,'b.');
    errorbar(r_arr.*1.01,profcb_all(im,:),profcb_all_err(im,:),'r.');
    xlim([4e-1,1.1e3])
    ylim([-1,3])
    if im==1
        ylabel('I [nW/m^2/sr]', 'fontsize',15);
    end
    if im==4
        h=legend({'PSF','HSC all','HSC all (field-averaged)'},...
            'Location','northeast');
        set(h,'fontsize',7)
    end
    subplot(3,4,im+8)
    semilogx(r_arr,psf_arr.*dat(im).profpsg(1),'k');hold on
    errorbar(r_arr.*0.99,dat(im).profpsg,dat(im).errps,'b.');
    errorbar(r_arr.*1.01,profps_all(im,:),profps_all_err(im,:),'r.');
    xlim([4e-1,1.1e3])
    ylim([-0.1,0.5])
    xlabel('arcsec', 'fontsize',15);
    if im==1
        ylabel('I [nW/m^2/sr]', 'fontsize',15);
    end
    if im==4
        h=legend({'PSF','HSC m<20','HSC m<20 (field-averaged)'},...
            'Location','northeast');
        set(h,'fontsize',7)
    end
end
%%
figure
setwinsize(gcf,1400,300)

for im=1:4
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    r_arr = dat(im).r_arr;
    
    subplot(1,4,im)

    loglog(r_arr, ecb_all(im,:),'k.','markersize',10);hold on
    loglog(r_arr, eps_all(im,:),'m.','markersize',10);hold on
    h=legend({'excess (all)','excess (m<20)'},...
        'Location','northeast');
    set(h,'fontsize',7)
    legend boxoff

    v = eps_all(im,:);
    e = profps_all_err(im,:);
    e(v-e<0) = v(v-e<0)-1e-8;
    loglog(r_arr, v,'m.','markersize',10);hold on
    errorbar(r_arr, v, e,'m.','markersize',10);
    v = -eps_all(im,:);
    e = profps_all_err(im,:);
    e(v-e<0) = v(v-e<0)-1e-8;    
    loglog(r_arr, v,'mo','markersize',5);hold on
    errorbar(r_arr, v, e,'mo','markersize',5);

    v = ecb_all(im,:);
    e = profcb_all_err(im,:);
    e(v-e<0) = v(v-e<0)-1e-8;
    loglog(r_arr, v,'k.','markersize',10);hold on
    errorbar(r_arr, v, e,'k.','markersize',10);
    v = -ecb_all(im,:);
    e = profcb_all_err(im,:);
    e(v-e<0) = v(v-e<0)-1e-8;    
    loglog(r_arr, v,'ko','markersize',5);hold on
    errorbar(r_arr, v, e,'ko','markersize',5);
    
    ylim([1e-4,1e3])
    xlim([4e-1,1.1e3])
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);
    xlabel('arcsec', 'fontsize',15);
    ylabel('I [nW/m^2/sr]', 'fontsize',15);

end
%%
%%%%%%%%%%%%%%%%%
%%
flight=40030;
mypaths=get_paths(flight);
inst=1;
masklim=false;
Niter = 2;
m_min_arr = 16:19;
m_max_arr = m_min_arr + 1;
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

ifield = 4;
dt=get_dark_times(flight,inst,ifield);


for im=1:4
    r_arr = stackdatfluc(im).stackdat.r_arr;
    excb = zeros([Niter,numel(r_arr)]);
    exps = zeros([Niter,numel(r_arr)]);
    excb100 = zeros([1,Niter]);
    exps100 = zeros([1,Niter]);
    profcbg = zeros([Niter,numel(r_arr)]);
    profcbs = zeros([Niter,numel(r_arr)]);
    profpsg = zeros([Niter,numel(r_arr)]);
    profpss = zeros([Niter,numel(r_arr)]);
    
for iter=1:6
    if masklim
        load(sprintf('%s/stackdatfluc_%s_masklim_iter%d',...
            loaddir,dt.name,iter+2),'stackdatfluc');        
    else
        load(sprintf('%s/stackdatfluc_%s_iter%d',...
            loaddir,dt.name,iter+2),'stackdatfluc');
    end
    
    excb(iter,:) = stackdatfluc(im).stackdat.excess.diffcb;
    exps(iter,:) = stackdatfluc(im).stackdat.excess.diffps;
    excb100(iter) = stackdatfluc(im).stackdat.excess.diffcb100;
    exps100(iter) = stackdatfluc(im).stackdat.excess.diffps100;
    profcbg(iter,:) = stackdatfluc(im).stackdat.all.profcbg;
    profcbs(iter,:) = stackdatfluc(im).stackdat.all.profcbs;
    profpsg(iter,:) = stackdatfluc(im).stackdat.all.profpsg;
    profpss(iter,:) = stackdatfluc(im).stackdat.all.profpss;
    
    counts = stackdatfluc(im).stackdat.all.counts;
    countg = stackdatfluc(im).stackdat.all.countg;
end
dat(im).r_arr = r_arr;
dat(im).counts = counts;
dat(im).countg = countg;
dat(im).excb_all = excb;
dat(im).exps_all = exps;
dat(im).excb100_all = excb100;
dat(im).exps100_all = exps100;
dat(im).excb = mean(excb);
dat(im).excb_err = std(excb);
dat(im).exps = mean(exps);
dat(im).exps_err = std(exps);
dat(im).excb100 = mean(excb100);
dat(im).excb100_err = std(excb100);
dat(im).exps100 = mean(exps100);
dat(im).exps_100err = std(exps100);
dat(im).profcbg = profcbg;
dat(im).profcbs = profcbs;
dat(im).profpsg = profpsg;
dat(im).profpss = profpss;

end
%%
figure
setwinsize(gcf,1400,360)
for im=1:4
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    r_arr = dat(im).r_arr;
    counts = dat(im).counts;
    countg = dat(im).countg;
    profcbg = dat(im).profcbg;
    profcbs = dat(im).profcbs;
    profpsg = dat(im).profpsg;
    profpss = dat(im).profpss;

    subplot(1,4,im)
    semilogx(r_arr.*1.01,profpss(1,:),'r.','markersize',10);hold on
    semilogx(r_arr.*0.98,profpsg(1,:),'b.','markersize',10);hold on
    semilogx(r_arr,profcbs(1,:),'r','markersize',10);hold on
    semilogx(r_arr,profcbg(1,:),'b','markersize',10);hold on
    h=legend({'PanSTARRS stars','PanSTARRS galaxies',...
        'PanSTARRS stars w/ fluctuations',...
        'PanSTARRS galaxies w/ fluctuations'},...
        'Location','northeast');
    set(h,'fontsize',7)
    legend boxoff

    for iter=1:6
        semilogx(r_arr,profcbg(iter,:),'b','markersize',10);hold on
        semilogx(r_arr,profcbs(iter,:),'r','markersize',10);hold on
    end
    xlim([1e1,1.1e3])
    ylim([-1,2])
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
      ' (',num2str(counts),' stars, ',num2str(countg), ' galaxies)'),'fontsize',10)
    xlabel('arcsec', 'fontsize',10)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)

end
if masklim
    suptitle(strcat(dt.name,' (mask to mag bin max)'));
else
    suptitle(strcat(dt.name,' (mask all PanSTARRS sources)'));
end
%%
figure
setwinsize(gcf,1400,360)
for im=1:4
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    r_arr = dat(im).r_arr;
    counts = dat(im).counts;
    countg = dat(im).countg;
    subplot(1,4,im)

    loglog(r_arr.*0.98, dat(im).excb,'r.','markersize',10);hold on
    loglog(r_arr.*1.02, dat(im).exps,'b.','markersize',10);hold on
    loglog(r_arr, dat(im).excb-dat(im).exps,'k.','markersize',10);hold on
    h=legend({'CIBER Excess','PanSTARRS Excess', 'Excess'},'Location','northeast');
    set(h,'fontsize',10)
    legend boxoff

    
    v = dat(im).excb;
    e = dat(im).excb_err;
    e(v-e<0) = v(v-e<0)-1e-8;
    loglog(r_arr.*0.98, v,'r.','markersize',10);hold on
    errorbar(r_arr.*0.98, v, e,'r.','markersize',10);
    v = -dat(im).excb;
    e = dat(im).excb_err;
    e(v-e<0) = v(v-e<0)-1e-8;    
    loglog(r_arr.*0.98, v,'ro','markersize',5);hold on
    errorbar(r_arr.*0.98, v, e,'ro','markersize',5);


    v = dat(im).exps;
    e = dat(im).exps_err;
    e(v-e<0) = v(v-e<0)-1e-8;
    loglog(r_arr.*1.02, v,'b.','markersize',10);hold on
    errorbar(r_arr.*1.02, v, e,'b.','markersize',10);
    v = -dat(im).exps;
    e = dat(im).exps_err;
    e(v-e<0) = v(v-e<0)-1e-8;    
    loglog(r_arr.*1.02, v,'bo','markersize',5);hold on
    errorbar(r_arr.*1.02, v, e,'bo','markersize',5);
    
    v = dat(im).excb - dat(im).exps;
    e = sqrt(dat(im).excb_err.^2 + dat(im).exps_err.^2);
    e(v-e<0) = v(v-e<0)-1e-8;
    loglog(r_arr, v,'k.','markersize',10);hold on
    errorbar(r_arr, v, e,'k.','markersize',10);
    v =- (dat(im).excb - dat(im).exps);
    e = sqrt(dat(im).excb_err.^2 + dat(im).exps_err.^2);
    e(v-e<0) = v(v-e<0)-1e-8;    
    loglog(r_arr, v,'ko','markersize',5);hold on
    errorbar(r_arr, v, e,'ko','markersize',5);

    xlim([4e-1,1.1e3])
    ylim([1e-4,1e4])
    title(strcat(num2str(m_min),'<m<',num2str(m_max),...
      ' (',num2str(counts),' stars, ',num2str(countg), ' galaxies)'),'fontsize',10)
    xlabel('arcsec', 'fontsize',10)
    ylabel('I [nW/m^2/sr]', 'fontsize',15)

end
if masklim
    suptitle(strcat(dt.name,' (mask to mag bin max)'));
else
    suptitle(strcat(dt.name,' (mask all PanSTARRS sources)'));
end
