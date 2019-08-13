flight=40030;
inst=1;
ifield=8;
hsc_idx=3;
Nsim=10;
rmin=nan;
sample_type = 'jack_random';
subpix=false;
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
name = HSC_fields_info(hsc_idx);
load(sprintf('%s/stackmapdathsc_%s%d',loaddir,name,0),...
    'stackmapdatsim');
dx = 1200;
verbose = 0;
cbmap0 = stackmapdatsim(ifield).all.cbmap;
psmap0 = stackmapdatsim(ifield).all.psmap;
mask_inst = stackmapdatsim(ifield).all.mask_inst_clip;
strmask = stackmapdatsim(ifield).all.strmask;
strnum = stackmapdatsim(ifield).all.strnum;

nbins = 25;
profile = radial_prof(ones(2*dx+1),ones(2*dx+1),dx+1,dx+1,1,nbins);
rbinedges = profile.binedges;
rbins = binedges2bins(rbinedges).*0.7;

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
%%
for isim=1:Nsim
    sigmap = sigmap_from_mz14(642,7);
    cbmap = cbmap0 + sigmap;
    psmap = cbmap0;
for im= 3:4

    m_min = im + 15;
    m_max = m_min + 1;
    
    stackdat.m_min = m_min;
    stackdat.m_max = m_max;
    
    srcdat = hsc_src_select(flight,inst,ifield,hsc_idx,m_min,m_max,...
        mask_inst,'sample_type',sample_type);
    
    clipmins = ones([2,nbins]).*-inf;
    clipmaxs = ones([2,nbins]).*inf;

    stackdat.r_arr = rbins;

    for isub=1:16
        [stampercbs,stamperpss,hitmaps,profcbs,profpss,profhits] = ...
            stackihl_ps0_hist_map(flight,inst,ifield,dx,cbmap,psmap,mask_inst,...
            strmask,strnum,1,verbose,rmin,clipmaxs,clipmins,...
            srcdat.sub(isub).xs_arr,srcdat.sub(isub).ys_arr,...
            srcdat.sub(isub).ms_arr,subpix);

        [stampercbg,stamperpsg,hitmapg,profcbg,profpsg,profhitg] = ...
            stackihl_ps0_hist_map(flight,inst,ifield,dx,cbmap,psmap,mask_inst,...
            strmask,strnum,1,verbose,rmin,clipmaxs,clipmins,...
            srcdat.sub(isub).xg_arr,srcdat.sub(isub).yg_arr,...
            srcdat.sub(isub).mg_arr,subpix);

        fprintf('stack isim %d, %d<m<%d, %s, isub %d\n',...
            isim,m_min,m_max,sample_type,isub);

        stackdat.sub(isub).counts = srcdat.sub(isub).Ns;
        stackdat.sub(isub).countg = srcdat.sub(isub).Ng;
        stackdat.sub(isub).profcbs = profcbs;
        stackdat.sub(isub).profcbg = profcbg;
        stackdat.sub(isub).profpss = profpss;
        stackdat.sub(isub).profpsg = profpsg;
        stackdat.sub(isub).profhits = profhits;
        stackdat.sub(isub).profhitg = profhitg;
    end
    sp100 = find(stackdat.r_arr>100);
  %% profile combining all subset
    profcbs = zeros(size(stackdat.r_arr));
    profcbg = zeros(size(stackdat.r_arr));
    profpss = zeros(size(stackdat.r_arr));
    profpsg = zeros(size(stackdat.r_arr));
    profhits = zeros(size(stackdat.r_arr));
    profhitg = zeros(size(stackdat.r_arr));
    counts = 0;
    countg = 0;
    for isub=1:16
        profcbs = profcbs + stackdat.sub(isub).profcbs.*stackdat.sub(isub).profhits;
        profcbg = profcbg + stackdat.sub(isub).profcbg.*stackdat.sub(isub).profhitg;
        profpss = profpss + stackdat.sub(isub).profpss.*stackdat.sub(isub).profhits;
        profpsg = profpsg + stackdat.sub(isub).profpsg.*stackdat.sub(isub).profhitg;
        profhits = profhits + stackdat.sub(isub).profhits;
        profhitg = profhitg + stackdat.sub(isub).profhitg;
        counts = counts + stackdat.sub(isub).counts;
        countg = countg + stackdat.sub(isub).countg;
    end
    stackdat.all.profcbs = profcbs./profhits;
    stackdat.all.profcbg = profcbg./profhitg;
    stackdat.all.profpss = profpss./profhits;
    stackdat.all.profpsg = profpsg./profhitg;
    stackdat.all.counts = counts;
    stackdat.all.countg = countg;
    
    stackdat.all.profcbs100 = sum(profcbs(sp100))./sum(profhits(sp100));
    stackdat.all.profcbg100 = sum(profcbg(sp100))./sum(profhitg(sp100));
    stackdat.all.profpss100 = sum(profpss(sp100))./sum(profhits(sp100));
    stackdat.all.profpsg100 = sum(profpsg(sp100))./sum(profhitg(sp100));

  %% profile of jackknife samples (leave one out)
    for isub=1:16
        jackcbs = profcbs - stackdat.sub(isub).profcbs.*stackdat.sub(isub).profhits;
        jackcbg = profcbg - stackdat.sub(isub).profcbg.*stackdat.sub(isub).profhitg;
        jackpss = profpss - stackdat.sub(isub).profpss.*stackdat.sub(isub).profhits;
        jackpsg = profpsg - stackdat.sub(isub).profpsg.*stackdat.sub(isub).profhitg;
        jackhits = profhits - stackdat.sub(isub).profhits;
        jackhitg = profhitg - stackdat.sub(isub).profhitg;
        stackdat.jack(isub).profcbs = jackcbs./jackhits;
        stackdat.jack(isub).profcbg = jackcbg./jackhitg;
        stackdat.jack(isub).profpss = jackpss./jackhits;
        stackdat.jack(isub).profpsg = jackpsg./jackhitg; 

        stackdat.jack(isub).profcbs100 = sum(jackcbs(sp100))./sum(jackhits);
        stackdat.jack(isub).profcbg100 = sum(jackcbg(sp100))./sum(jackhitg);
        stackdat.jack(isub).profpss100 = sum(jackpss(sp100))./sum(jackhits);
        stackdat.jack(isub).profpsg100 = sum(jackpsg(sp100))./sum(jackhitg); 
        
    end
  %% error bar with jackknife
    errcbs = zeros(size(stackdat.r_arr));
    errcbg = zeros(size(stackdat.r_arr));
    errpss = zeros(size(stackdat.r_arr));
    errpsg = zeros(size(stackdat.r_arr));

    errcbs100 = 0;
    errcbg100 = 0;
    errpss100 = 0;
    errpsg100 = 0;

    for isub=1:16
    errcbs = errcbs + (stackdat.jack(isub).profcbs - stackdat.all.profcbs).^2;
    errcbg = errcbg + (stackdat.jack(isub).profcbg - stackdat.all.profcbg).^2;
    errpss = errpss + (stackdat.jack(isub).profpss - stackdat.all.profpss).^2;
    errpsg = errpsg + (stackdat.jack(isub).profpsg - stackdat.all.profpsg).^2;

    errcbs100 = errcbs100 + ...
        (stackdat.jack(isub).profcbs100 - stackdat.all.profcbs100).^2;
    errcbg100 = errcbg100 + ...
        (stackdat.jack(isub).profcbg100 - stackdat.all.profcbg100).^2;
    errpss100 = errpss100 + ...
        (stackdat.jack(isub).profpss100 - stackdat.all.profpss100).^2;
    errpsg100 = errpsg100 + ...
        (stackdat.jack(isub).profpsg100 - stackdat.all.profpsg100).^2;
    end
    stackdat.errjack.profcbs = sqrt(errcbs.*(15/16));
    stackdat.errjack.profcbg = sqrt(errcbg.*(15/16));
    stackdat.errjack.profpss = sqrt(errpss.*(15/16));
    stackdat.errjack.profpsg = sqrt(errpsg.*(15/16));
    stackdat.errjack.profcbs100 = sqrt(errcbs100.*(15/16));
    stackdat.errjack.profcbg100 = sqrt(errcbg100.*(15/16));
    stackdat.errjack.profpss100 = sqrt(errpss100.*(15/16));
    stackdat.errjack.profpsg100 = sqrt(errpsg100.*(15/16));

    
%%    
    dat(im).isim(isim).stackdat = stackdat;
    clear stackdat
end
end
dathsc(hsc_idx).dat = dat;

%%
hsc_idx=2;
hsc_name = HSC_fields_info(hsc_idx);
r_arr = dathsc(1).dat(3).isim(1).stackdat.r_arr;

figure
setwinsize(gcf,800,600)
for im=3:4
    ecb_arr = zeros([Nsim,numel(r_arr)]);
    eps_arr = zeros([Nsim,numel(r_arr)]);
    m_min = dathsc(hsc_idx).dat(im).isim(1).stackdat.m_min;
    m_max = dathsc(hsc_idx).dat(im).isim(1).stackdat.m_max;    
    for isim=1:Nsim
        m_min = dathsc(hsc_idx).dat(im).isim(isim).stackdat.m_min;
        m_max = dathsc(hsc_idx).dat(im).isim(isim).stackdat.m_max;

        profcbs = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profcbs;
        profcbg = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profcbg;
        ecb_arr(isim,:) = profcbg - profcbs;
        profpss = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profpss;
        profpsg = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profpsg;
        eps_arr(isim,:) = profpsg - profpss;
    end
    subplot(2,2,im-2)
    semilogx(r_arr,nanmean(ecb_arr),'r.-');hold on
    semilogx(r_arr,nanmean(eps_arr),'b.-');hold on
    legend({'src + fluc','src'});
    errorbar(r_arr, nanmean(ecb_arr), nanstd(ecb_arr),'r.');
    xlim([1e1,1.1e3])
    ylim([-0.4,0.4])
    hline(0,'k--')
    xlabel('arcsec', 'fontsize',10)
    ylabel('gal stack -star stack [nW/m^2/sr]', 'fontsize',10)
    title(sprintf('%d<m<%d one field',m_min,m_max));
end

for im=3:4
    ecb_arr = zeros([Nsim*2,numel(r_arr)]);
    eps_arr = zeros([Nsim*2,numel(r_arr)]);
    m_min = dathsc(hsc_idx).dat(im).isim(1).stackdat.m_min;
    m_max = dathsc(hsc_idx).dat(im).isim(1).stackdat.m_max;
    count=0;
    for hsc_idx=1:3
        for isim=1:Nsim
            count=count+1;
            m_min = dathsc(hsc_idx).dat(im).isim(isim).stackdat.m_min;
            m_max = dathsc(hsc_idx).dat(im).isim(isim).stackdat.m_max;

            profcbs = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profcbs;
            profcbg = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profcbg;
            ecb_arr(count,:) = profcbg - profcbs;
            profpss = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profpss;
            profpsg = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profpsg;
            eps_arr(count,:) = profpsg - profpss;
        end
    end
    subplot(2,2,im)
    semilogx(r_arr.*1.01,nanmean(ecb_arr),'r.-');hold on
    semilogx(r_arr.*0.99,nanmean(eps_arr),'b.-');
    legend({'src + fluc','src'});
    errorbar(r_arr.*1.01, nanmean(ecb_arr), nanstd(ecb_arr),'r.');
    errorbar(r_arr.*0.99, nanmean(eps_arr), nanstd(eps_arr),'b.');
    xlim([1e1,1.1e3])
    ylim([-0.4,0.4])
    hline(0,'k--')
    xlabel('arcsec', 'fontsize',10)
    ylabel('gal stack -star stack [nW/m^2/sr]', 'fontsize',10)
    title(sprintf('%d<m<%d all fields',m_min,m_max));
end
print(sprintf('/Users/ytcheng/Desktop/ciber/fluc_test'),'-dpng');

%% plot individual profile
hsc_idx=2;
hsc_name = HSC_fields_info(hsc_idx);
r_arr = dathsc(1).dat(3).isim(1).stackdat.r_arr;
for isim=1:Nsim
    figure
    setwinsize(gcf,800,600)
    for im=3:4
        m_min = dathsc(hsc_idx).dat(im).isim(isim).stackdat.m_min;
        m_max = dathsc(hsc_idx).dat(im).isim(isim).stackdat.m_max;

        profcbs = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profcbs;
        profcbg = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profcbg;
        errcbs = dathsc(hsc_idx).dat(im).isim(isim).stackdat.errjack.profcbs;
        errcbg = dathsc(hsc_idx).dat(im).isim(isim).stackdat.errjack.profcbg;

        profpss = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profpss;
        profpsg = dathsc(hsc_idx).dat(im).isim(isim).stackdat.all.profpsg;
        errpss = dathsc(hsc_idx).dat(im).isim(isim).stackdat.errjack.profpss;
        errpsg = dathsc(hsc_idx).dat(im).isim(isim).stackdat.errjack.profpsg;

        subplot(2,2,im-2)
        semilogx(r_arr.*1.01, profcbs, 'r');hold on
        semilogx(r_arr.*0.99, profcbg, 'b');
        legend({'stars','galaxies'});
        hline(0,'k--');
        errorbar(r_arr.*1.01, profcbs, errcbs,'r.');
        errorbar(r_arr.*0.99, profcbg, errcbg,'b.');
        title(sprintf('%d<m<%d src + fluc',m_min,m_max));
        xlim([1e1,1.1e3])
        ylim([-0.8,1])

        subplot(2,2,im)
        semilogx(r_arr.*1.01, profpss, 'r');hold on
        semilogx(r_arr.*0.99, profpsg, 'b');
        legend({'stars','galaxies'});
        hline(0,'k--');
        errorbar(r_arr.*1.01, profpss, errpss,'r.');
        errorbar(r_arr.*0.99, profpsg, errpsg,'b.');
        title(sprintf('%d<m<%d fluc',m_min,m_max));
        xlim([1e1,1.1e3])
        ylim([-0.08,0.1])
    end
    print(sprintf('/Users/ytcheng/Desktop/ciber/%s_fluc_test%d',...
    hsc_name,isim),'-dpng');
end