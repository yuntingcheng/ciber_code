function stack_preprocess_sim(flight,inst,f_ihl,rvir)

mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'maskdat'),'maskdat');

savedir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackmapdat',savedir),'stackmapdat');

srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');

for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);

    cbmap = stackmapdat(ifield).cbmap;
    cbmask = stackmapdat(ifield).strmask.*stackmapdat(ifield).mask_inst_clip; 
    sign = std(cbmap(cbmask==1));

    mask_inst = ones(720);
    strmask = maskdat.mask(ifield).strmask_sim;
    totmask = mask_inst.*strmask;
    
    map1 = fits_read(strcat(srcmapdir,dt.name,'_srcmap_sim1_all.fits'));
    map2 = fits_read(strcat(srcmapdir,dt.name,'_srcmap_sim2_all.fits'));
    map3 = fits_read(strcat(srcmapdir,dt.name,'_srcmap_sim3_all.fits'));
    
    ihlmap1 = fits_read(strcat(srcmapdir,...
        'unisphere_ihlmap_sim1_all',num2str(rvir),'.fits'));
    ihlmap2 = fits_read(strcat(srcmapdir,...
        'unisphere_ihlmap_sim2_all',num2str(rvir),'.fits'));

    psmap_raw = map1;

    %%%%%%%%%%%%%%%%%%%% all sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cbmap_raw = map1 + map2 + map3 + ihlmap1.*f_ihl + ihlmap2.*f_ihl;
    %%% sigma clip and mean/grad sub %%%
%     Q1 = quantile(cbmap_raw(find(totmask)),0.25);
%     Q3 = quantile(cbmap_raw(find(totmask)),0.75);
%     IQR = Q3-Q1;
%     clipmin = Q1 - 30*IQR;
%     clipmax = Q3 + 30*IQR;
%     sigmask = totmask;
%     sigmask((cbmap_raw>clipmax) | (cbmap_raw<clipmin)) = 0;
    sigmask = totmask;
    
    cbmean = mean(cbmap_raw(find(sigmask)));
    psmean = mean(psmap_raw(find(sigmask)));
    
    p = polyfitweighted2(1:720,1:720,cbmap_raw,2,sigmask);
    polymapcb = polyval2(p,1:720,1:720);
    cbmap = cbmap_raw - polymapcb;
    cb_bk = mean(cbmap(find(sigmask)));
    cbmap = cbmap - cb_bk;

    p = polyfitweighted2(1:720,1:720,psmap_raw,2,sigmask);
    polymapps = polyval2(p,1:720,1:720);
    psmap = psmap_raw - polymapps;
    ps_bk = mean(psmap(find(sigmask)));
    psmap = psmap - ps_bk;

    sig_sp = find((totmask-sigmask)==1);
    mask_inst_clip = mask_inst;
    mask_inst_clip(sig_sp)=0;

    %%% write the data %%%
    stackmapdatsim(ifield).all.cbmap = cbmap;
    stackmapdatsim(ifield).all.psmap = psmap;
    stackmapdatsim(ifield).all.mask_inst_clip = mask_inst_clip;
    stackmapdatsim(ifield).all.strmask = maskdat.mask(ifield).strmask_sim;
    stackmapdatsim(ifield).all.strnum = maskdat.mask(ifield).strnum_sim;
    stackmapdatsim(ifield).all.cbmean = cbmean;
    stackmapdatsim(ifield).all.psmean = psmean;
    stackmapdatsim(ifield).all.sign = sign;

    %%%%%%%%%%%%%%%%%%%% partial sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cbmap_raw = map1 + map2 + ihlmap1.*f_ihl + ihlmap2.*f_ihl;
    
    %%% sigma clip and mean/grad sub %%%
%     Q1 = quantile(cbmap_raw(find(totmask)),0.25);
%     Q3 = quantile(cbmap_raw(find(totmask)),0.75);
%     IQR = Q3-Q1;
%     clipmin = Q1 - 30*IQR;
%     clipmax = Q3 + 30*IQR;
%     sigmask = totmask;
%     sigmask((cbmap_raw>clipmax) | (cbmap_raw<clipmin)) = 0;
    sigmask = totmask;
    
    cbmean = mean(cbmap_raw(find(sigmask)));
    psmean = mean(psmap_raw(find(sigmask)));
    
    p = polyfitweighted2(1:720,1:720,cbmap_raw,2,sigmask);
    polymapcb = polyval2(p,1:720,1:720);
    cbmap = cbmap_raw - polymapcb;
    cb_bk = mean(cbmap(find(sigmask)));
    cbmap = cbmap - cb_bk;

    p = polyfitweighted2(1:720,1:720,psmap_raw,2,sigmask);
    polymapps = polyval2(p,1:720,1:720);
    psmap = psmap_raw - polymapps;
    ps_bk = mean(psmap(find(sigmask)));
    psmap = psmap - ps_bk;

    sig_sp = find((totmask-sigmask)==1);
    mask_inst_clip = mask_inst;
    mask_inst_clip(sig_sp)=0;

    %%% write the data %%%
    stackmapdatsim(ifield).sub.cbmap = cbmap;
    stackmapdatsim(ifield).sub.psmap = psmap;
    stackmapdatsim(ifield).sub.mask_inst_clip = mask_inst_clip;
    stackmapdatsim(ifield).sub.strmask = maskdat.mask(ifield).strmask_sim;
    stackmapdatsim(ifield).sub.strnum = maskdat.mask(ifield).strnum_sim;
    stackmapdatsim(ifield).sub.cbmean = cbmean;
    stackmapdatsim(ifield).sub.psmean = psmean;
    stackmapdatsim(ifield).sub.sign = sign;
    stackmapdatsim(ifield).sub.m_max = 22;

end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
if rvir==1
save(sprintf('%s/stackmapdatsim%d',savedir,f_ihl*100),'stackmapdatsim');
else
save(sprintf('%s/stackmapdatsim%d_rv%d',savedir,f_ihl*100,rvir),'stackmapdatsim');
end
return