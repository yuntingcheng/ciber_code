function stack_preprocess_sim(flight,inst)

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
    
    psmap_raw = map1;

    %%%%%%%%%%%%%%%%%%%% all sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cbmap_raw = map1 + map2 + map3 + normrnd(0, sign, size(map1));
    
    %%% sigma clip and mean/grad sub %%%
    sigmask1 = sigclip_mask(cbmap_raw,totmask,3,5);
    sigmask1 = sigclip_mask(psmap_raw,sigmask1,3,5);
    sm = fillpadsmooth(cbmap_raw,sigmask1,2);
    sigmask2 = sigclip_mask(sm,sigmask1,3,5);
    sm = fillpadsmooth(psmap_raw,sigmask2,2);
    sigmask = sigclip_mask(sm,sigmask2,3,5);
    
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
    cbmap_raw = map1 + map2 + normrnd(0, sign, size(map1));
    
    %%% sigma clip and mean/grad sub %%%
    sigmask1 = sigclip_mask(cbmap_raw,totmask,3,5);
    sigmask1 = sigclip_mask(psmap_raw,sigmask1,3,5);
    sm = fillpadsmooth(cbmap_raw,sigmask1,2);
    sigmask2 = sigclip_mask(sm,sigmask1,3,5);
    sm = fillpadsmooth(psmap_raw,sigmask2,2);
    sigmask = sigclip_mask(sm,sigmask2,3,5);
    
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
save(sprintf('%s/stackmapdatsim',savedir),'stackmapdatsim');

return