function run_ihl_stack_map_fluc(flight,inst,ifield,iter_arr,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addRequired('iter_arr',@isnumeric);
  p.addOptional('masklim',false,@islogical);
  p.addOptional('sample_type','jack_random',@ischar);

  p.parse(flight,inst,ifield,iter_arr,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  iter_arr   = p.Results.iter_arr;
  masklim = p.Results.masklim;
  sample_type=p.Results.sample_type;
  clear p varargin;
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
if masklim
    load(sprintf('%s/stackdat_%s_masklim',...
        loaddir,dt.name),'stackdatall');        
else
    load(sprintf('%s/stackdat_%s',...
        loaddir,dt.name),'stackdatall');
end

dx = 1200;
verbose = false;
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;
m_min_arr = m_min_arr(10:13);
m_max_arr = m_max_arr(10:13);
nsim = 30;
if inst == 1
    stackband = 'I';
else
    stackband = 'H';
end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

nbins = 25;
profile = radial_prof(ones(2*dx+1),ones(2*dx+1),dx+1,dx+1,1,nbins);
rbinedges = profile.binedges;
rbins = binedges2bins(rbinedges).*0.7;
sp100 = find(rbins>100);
%%
for im= 1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    if masklim
        mask_inst = stackmapdat(ifield).m_max(m_max).mask_inst_clip;
        strmask = stackmapdat(ifield).m_max(m_max).strmask;
        strnum = stackmapdat(ifield).m_max(m_max).strnum;
    else
        mask_inst = stackmapdat(ifield).mask_inst_clip;
        strmask = stackmapdat(ifield).strmask;
        strnum = stackmapdat(ifield).strnum;    
    end
    srcdat = ps_src_select(flight,inst,ifield,m_min,m_max,mask_inst,stackband,...
    'sample_type',sample_type);

    [clipmaxs, clipmins, ~]=...
    stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
    mask_inst,strnum,1000,verbose,stackband,[],nan);
    
    tempdat(im).mask_inst = mask_inst;
    tempdat(im).strmask = strmask;
    tempdat(im).strnum = strnum;
    tempdat(im).srcdat = srcdat;
    tempdat(im).clipmaxs = clipmaxs;
    tempdat(im).clipmins = clipmins;
    tempdat(im).normcb = stackdatall(im).stackdat.norm.normcb;
    tempdat(im).normps = stackdatall(im).stackdat.norm.normps;
    
end
%%
for iter=iter_arr
cbmap1 = psmap + sigmap_from_mz14(1024,7);
for im= 1:numel(m_min_arr)
    stackdat.r_arr = rbins;
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    stackdat.m_min = m_min;
    stackdat.m_max = m_max;
    mask_inst = tempdat(im).mask_inst;
    strmask = tempdat(im).strmask;
    strnum = tempdat(im).strnum;
    srcdat = tempdat(im).srcdat;
    clipmaxs = tempdat(im).clipmaxs;
    clipmins = tempdat(im).clipmins;

    %%% background stack %%%
    r_arr = rbins;
    profcbs_arr = zeros([nsim,numel(r_arr)]);
    profcbg_arr = zeros([nsim,numel(r_arr)]);
    profpss_arr = zeros([nsim,numel(r_arr)]);
    profpsg_arr = zeros([nsim,numel(r_arr)]);
    profcbs100 = zeros([1,nsim]);
    profcbg100 = zeros([1,nsim]);
    profpss100 = zeros([1,nsim]);
    profpsg100 = zeros([1,nsim]);
    for isim=1:nsim
        [profcb,profps,hitmap]=stackihl_ps0_hist_map_bk...
            (dx,cbmap1,psmap,mask_inst,strmask,[srcdat.Ns,srcdat.Ng],...
            verbose,false);
        profcbs_arr(isim,:) = profcb(1,:);
        profcbg_arr(isim,:) = profcb(2,:);
        profpss_arr(isim,:) = profps(1,:);
        profpsg_arr(isim,:) = profps(2,:);
        profcbs100(isim) = sum(profcb(1,sp100).*hitmap(1,sp100))./...
            sum(hitmap(1,sp100));
        profcbg100(isim) = sum(profcb(2,sp100).*hitmap(2,sp100))./...
            sum(hitmap(2,sp100));
        profpss100(isim) = sum(profps(1,sp100).*hitmap(1,sp100))./...
            sum(hitmap(1,sp100));
        profpsg100(isim) = sum(profps(2,sp100).*hitmap(2,sp100))./...
            sum(hitmap(2,sp100));
        fprintf('stack %s, iter %d, %d<m<%d, %d stars %d gals, isim %d\n',...
            dt.name,iter,m_min,m_max,srcdat.Ns,srcdat.Ng,isim);
    end
    profcbs_err = nanstd(profcbs_arr);
    profcbg_err = nanstd(profcbg_arr);
    profpss_err = nanstd(profpss_arr);
    profpsg_err = nanstd(profpsg_arr);
    profcbs_arr = nanmean(profcbs_arr);
    profcbg_arr = nanmean(profcbg_arr);
    profpss_arr = nanmean(profpss_arr);
    profpsg_arr = nanmean(profpsg_arr);
    sp = find(profcbs_arr==profcbs_arr);
    profcbs_arr = spline(r_arr(sp),profcbs_arr(sp),r_arr);
    sp = find(profcbg_arr==profcbg_arr);
    profcbg_arr = spline(r_arr(sp),profcbg_arr(sp),r_arr);
    sp = find(profpss_arr==profpss_arr);
    profpss_arr = spline(r_arr(sp),profpss_arr(sp),r_arr);
    sp = find(profpsg_arr==profpsg_arr);
    profpsg_arr = spline(r_arr(sp),profpsg_arr(sp),r_arr);
    sp = find(profcbs_err==profcbs_err);
    profcbs_err = spline(r_arr(sp),profcbs_err(sp),r_arr);
    sp = find(profcbg_err==profcbg_err);
    profcbg_err = spline(r_arr(sp),profcbg_err(sp),r_arr);
    sp = find(profpss_err==profpss_err);
    profpss_err = spline(r_arr(sp),profpss_err(sp),r_arr);
    sp = find(profpsg_err==profpsg_err);
    profpsg_err = spline(r_arr(sp),profpsg_err(sp),r_arr);  
    profcbs_100err = nanstd(profcbs100);
    profcbg_100err = nanstd(profcbg100);
    profpss_100err = nanstd(profpss100);
    profpsg_100err = nanstd(profpsg100);
    profcbs_100 = nanmean(profcbs100);
    profcbg_100 = nanmean(profcbg100);
    profpss_100 = nanmean(profpss100);
    profpsg_100 = nanmean(profpsg100);

    stackdat.bk.profcbs = profcbs_arr;
    stackdat.bk.profcbs_err= profcbs_err;
    stackdat.bk.profcbg = profcbg_arr;
    stackdat.bk.profcbg_err= profcbg_err;
    stackdat.bk.profpss = profpss_arr;
    stackdat.bk.profpss_err= profpss_err;
    stackdat.bk.profpsg = profpsg_arr;
    stackdat.bk.profpsg_err= profpsg_err;
    stackdat.bk.profcbs100 = profcbs_100;
    stackdat.bk.profcbs_err100= profcbs_100err;
    stackdat.bk.profcbg100 = profcbg_100;
    stackdat.bk.profcbg_err100= profcbg_100err;
    stackdat.bk.profpss100 = profpss_100;
    stackdat.bk.profpss_err100= profpss_100err;
    stackdat.bk.profpsg100 = profpsg_100;
    stackdat.bk.profpsg_err100= profpsg_100err;
    %%%%%%%%%
    
    for isub=1:16
        [~,~,~,profcbs,profpss,profhits] = ...
            stackihl_ps0_hist_map(flight,inst,ifield,dx,cbmap1,psmap,mask_inst,...
            strmask,strnum,1,verbose,nan,clipmaxs,clipmins,...
            srcdat.sub(isub).xs_arr,srcdat.sub(isub).ys_arr,...
            srcdat.sub(isub).ms_arr,false);

        [~,~,~,profcbg,profpsg,profhitg] = ...
            stackihl_ps0_hist_map(flight,inst,ifield,dx,cbmap1,psmap,mask_inst,...
            strmask,strnum,1,verbose,nan,clipmaxs,clipmins,...
            srcdat.sub(isub).xg_arr,srcdat.sub(isub).yg_arr,...
            srcdat.sub(isub).mg_arr,false);

        fprintf('stack %s, iter %d, %d<m<%d, isub %d\n',...
            dt.name,iter, m_min,m_max,isub);

        profcbs(profhits==0) = 0;
        profpss(profhits==0) = 0;
        profcbg(profhitg==0) = 0;
        profpsg(profhitg==0) = 0;

        stackdat.sub(isub).counts = srcdat.sub(isub).Ns;
        stackdat.sub(isub).countg = srcdat.sub(isub).Ng;
        stackdat.sub(isub).profcbs = profcbs;
        stackdat.sub(isub).profcbg = profcbg;
        stackdat.sub(isub).profpss = profpss;
        stackdat.sub(isub).profpsg = profpsg;
        stackdat.sub(isub).profhits = profhits;
        stackdat.sub(isub).profhitg = profhitg;
    end
   
    %%% profile combining all subset
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

    %%% profile of jackknife samples (leave one out)
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
    %%% error bar with jackknife
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

        
    %%% interpolate missing points
    p = stackdat.all.profcbs;
    stackdat.all.profcbs = spline(r_arr(find(p==p)),p(find(p==p)),r_arr);
    p = stackdat.all.profcbg;
    stackdat.all.profcbg = spline(r_arr(find(p==p)),p(find(p==p)),r_arr);
    p = stackdat.all.profpss;
    stackdat.all.profpss = spline(r_arr(find(p==p)),p(find(p==p)),r_arr);
    p = stackdat.all.profpsg;
    stackdat.all.profpsg = spline(r_arr(find(p==p)),p(find(p==p)),r_arr);

    p = stackdat.errjack.profcbs;
    stackdat.errjack.profcbs = spline(r_arr(find(p==p)),p(find(p==p)),r_arr);
    p = stackdat.errjack.profcbg;
    stackdat.errjack.profcbg = spline(r_arr(find(p==p)),p(find(p==p)),r_arr);
    p = stackdat.errjack.profpss;
    stackdat.errjack.profpss = spline(r_arr(find(p==p)),p(find(p==p)),r_arr);
    p = stackdat.errjack.profpsg;
    stackdat.errjack.profpsg = spline(r_arr(find(p==p)),p(find(p==p)),r_arr);
    
    %%% get normalized profile
    profcbg = stackdat.all.profcbg - stackdat.bk.profcbg;
    profcbg100 = stackdat.all.profcbg100 - stackdat.bk.profcbg100;
    profcbg_err = sqrt(stackdat.errjack.profcbg.^2 + ...
        stackdat.bk.profcbg_err.^2);
    profcbg_err100 = sqrt(stackdat.errjack.profcbg100.^2 + ...
        stackdat.bk.profcbg_err100.^2);
    
    profcbs = stackdat.all.profcbs - stackdat.bk.profcbs;
    profcbs100 = stackdat.all.profcbs100 - stackdat.bk.profcbs100;
    profcbs_err = sqrt(stackdat.errjack.profcbs.^2 + ...
        stackdat.bk.profcbs_err.^2);
    profcbs_err100 = sqrt(stackdat.errjack.profcbs100.^2 + ...
        stackdat.bk.profcbs_err100.^2);
    norm = tempdat(im).normcb;
    profcbs = profcbs.*norm;
    profcbs100 = profcbs100.*norm;
    profcbs_err = profcbs_err.*norm;
    profcbs_err100 = profcbs_err100.*norm;
    
    stackdat.norm.profcbg = profcbg;
    stackdat.norm.profcbg100 = profcbg100;
    stackdat.norm.profcbg_err = profcbg_err;
    stackdat.norm.profcbg_err100 = profcbg_err100;
    stackdat.norm.profcbs = profcbs;
    stackdat.norm.profcbs100 = profcbs100;
    stackdat.norm.profcbs_err = profcbs_err;
    stackdat.norm.profcbs_err100 = profcbs_err100;

    profpsg = stackdat.all.profpsg - stackdat.bk.profpsg;
    profpsg100 = stackdat.all.profpsg100 - stackdat.bk.profpsg100;
    profpsg_err = sqrt(stackdat.errjack.profpsg.^2 + ...
        stackdat.bk.profpsg_err.^2);
    profpsg_err100 = sqrt(stackdat.errjack.profpsg100.^2 + ...
        stackdat.bk.profpsg_err100.^2);
    
    profpss = stackdat.all.profpss - stackdat.bk.profpss;
    profpss100 = stackdat.all.profpss100 - stackdat.bk.profpss100;
    profpss_err = sqrt(stackdat.errjack.profpss.^2 + ...
        stackdat.bk.profpss_err.^2);
    profpss_err100 = sqrt(stackdat.errjack.profpss100.^2 + ...
        stackdat.bk.profpss_err100.^2);
    norm = tempdat(im).normps;
    profpss = profpss.*norm;
    profpss100 = profpss100.*norm;
    profpss_err = profpss_err.*norm;
    profpss_err100 = profpss_err100.*norm;

    stackdat.norm.profpsg = profpsg;
    stackdat.norm.profpsg100 = profpsg100;
    stackdat.norm.profpsg_err = profpsg_err;
    stackdat.norm.profpsg_err100 = profpsg_err100;
    stackdat.norm.profpss = profpss;
    stackdat.norm.profpss100 = profpss100;
    stackdat.norm.profpss_err = profpss_err;
    stackdat.norm.profpss_err100 = profpss_err100;

    %%% get excess profile
    diffcb = stackdat.norm.profcbg - stackdat.norm.profcbs;
    diffcb100 = stackdat.norm.profcbg100 - stackdat.norm.profcbs100;
    diffcb_err = sqrt(stackdat.norm.profcbg_err.^2 + ...
        stackdat.norm.profcbs_err.^2);
    diffcb_err100 = sqrt(stackdat.norm.profcbg_err100.^2 + ...
        stackdat.norm.profcbs_err100.^2);
    diffps = stackdat.norm.profpsg - stackdat.norm.profpss;
    diffps100 = stackdat.norm.profpsg100 - stackdat.norm.profpss100;
    diffps_err = sqrt(stackdat.norm.profpsg_err.^2 + ...
        stackdat.norm.profpss_err.^2);
    diffps_err100 = sqrt(stackdat.norm.profpsg_err100.^2 + ...
        stackdat.norm.profpss_err100.^2);
    diff = diffcb - diffps;
    diff100 = diffcb100 - diffps100;
    diff_err = sqrt(diffcb_err.^2 + diffps_err.^2);
    diff_err100 = sqrt(diffcb_err100.^2 + diffps_err100.^2);
    diffcb(1) = 0; diffps(1) = 0; diff(1) = 0;
    
    stackdat.excess.diffcb = diffcb;
    stackdat.excess.diffcb100 = diffcb100;
    stackdat.excess.diffcb_err = diffcb_err;
    stackdat.excess.diffcb_err100 = diffcb_err100;
    stackdat.excess.diffps = diffps;
    stackdat.excess.diffps100 = diffps100;
    stackdat.excess.diffps_err = diffps_err;
    stackdat.excess.diffps_err100 = diffps_err100;
    stackdat.excess.diff = diff;
    stackdat.excess.diff100 = diff100;
    stackdat.excess.diff_err = diff_err;
    stackdat.excess.diff_err100 = diff_err100;
    
    stackdatfluc(im).stackdat = stackdat;
    clear stackdat
    
    if masklim
        save(sprintf('%s/stackdatfluc_%s_masklim_iter%d',...
            savedir,dt.name,iter),'stackdatfluc');        
    else
        save(sprintf('%s/stackdatfluc_%s_iter%d',...
            savedir,dt.name,iter),'stackdatfluc');
    end

end
    clear stackdatfluc
end

return
