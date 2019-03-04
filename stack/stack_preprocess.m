function stack_preprocess(flight,inst)

mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'maskdat'),'maskdat');

m_min_arr = [0,8:22];
m_max_arr = [8:23];

savedir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackmapdat',savedir),'stackmapdat');

cal = get_cal_apf2nWpm2ps(inst);
for ifield=4:8

    dt=get_dark_times(flight,inst,ifield);
    cbmap_raw = stackmapdat(ifield).map * cal(ifield).apf2nWpm2ps;

    %%%%% PanSTARRS srcmap %%%%%%%
    srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
        num2str(inst),'/');
    psmap_raw = fits_read(strcat(srcmapdir,dt.name,'_srcmap_ps1_all.fits'));
    
    %%% masks %%%
    mask_inst = stackmapdat(ifield).mask;
    strmask = maskdat.mask(ifield).strmask_stack;
%     if ifield== 5 
%         strmask = maskdat.mask(ifield).strmask_stack_aggressive;
%     end
    totmask = mask_inst.*strmask;

    %%% sigma clip and mean/grad sub %%%
    sigmask1 = sigclip_mask(cbmap_raw,totmask,3,5);
    sigmask1 = sigclip_mask(psmap_raw,sigmask1,3,5);
    sm = fillpadsmooth(cbmap_raw,sigmask1,2);
    sigmask2 = sigclip_mask(sm,sigmask1,3,5);
    sm = fillpadsmooth(psmap_raw,sigmask2,2);
    sigmask = sigclip_mask(sm,sigmask2,3,5);
    
    cbmean = mean(cbmap_raw(find(sigmask)));
    psmean = mean(psmap_raw(find(sigmask)));
    
    p = polyfitweighted2(1:1024,1:1024,cbmap_raw,2,sigmask);
    polymapcb = polyval2(p,1:1024,1:1024);
    cbmap = cbmap_raw - polymapcb;
    cb_bk = mean(cbmap(find(sigmask)));
    cbmap = cbmap - cb_bk;
    
    p = polyfitweighted2(1:1024,1:1024,psmap_raw,2,sigmask);
    polymapps = polyval2(p,1:1024,1:1024);
    psmap = psmap_raw - polymapps;
    ps_bk = mean(psmap(find(sigmask)));
    psmap = psmap - ps_bk;

    sig_sp = find((totmask-sigmask)==1);
    mask_inst_clip = mask_inst;
    mask_inst_clip(sig_sp)=0;

    %%% get smoothed FF err %%%
    sm = fillpadsmooth(cbmap,mask_inst_clip.*strmask,50);
    
    %%% write the data %%%
    stackmapdat(ifield).cbmap = cbmap;
    stackmapdat(ifield).psmap = psmap;
    stackmapdat(ifield).mask_inst_clip = mask_inst_clip;
    stackmapdat(ifield).strmask = maskdat.mask(ifield).strmask_stack;
    stackmapdat(ifield).strnum = maskdat.mask(ifield).strnum_stack;
%     if ifield== 5 % 
%         stackmapdat(ifield).strmask = maskdat.mask(ifield).strmask_stack_aggressive;
%         stackmapdat(ifield).strnum = maskdat.mask(ifield).strnum_stack_aggressive;
%     end
    stackmapdat(ifield).psmap_FFerr = psmap + sm;
    stackmapdat(ifield).m_min_arr = m_min_arr;
    stackmapdat(ifield).m_max_arr = m_max_arr;
    stackmapdat(ifield).cbmean = cbmean;
    stackmapdat(ifield).psmean = psmean;
end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/stackmapdat',savedir),'stackmapdat');

return