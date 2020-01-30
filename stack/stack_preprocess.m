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
    psmap_raw = fits_read(strcat(srcmapdir,dt.name,'_srcmap_ps_all.fits'));
    
    %%% masks %%%
    mask_inst = stackmapdat(ifield).mask;
    strmask = maskdat.mask(ifield).m_max(20).strmask_stack;    
    totmask = mask_inst.*strmask;
    N1 = sum(totmask(:));%%%
    
    %%% sigma clip and mean/grad sub %%%
    Q1 = quantile(cbmap_raw(find(totmask)),0.25);
    Q3 = quantile(cbmap_raw(find(totmask)),0.75);
    IQR = Q3-Q1;
    clipmin = Q1 - 3*IQR;
    clipmax = Q3 + 3*IQR;
    sigmask = totmask;
    sigmask((cbmap_raw>clipmax) | (cbmap_raw<clipmin)) = 0;
    
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

    %%% clip the residual point sources in the map %%%
    smcb = fillpadsmooth(cbmap.*mask_inst_clip.*strmask,mask_inst.*strmask,1);
    Q1 = quantile(cbmap(find(smcb.*mask_inst_clip.*strmask)),0.25);
    Q3 = quantile(cbmap(find(smcb.*mask_inst_clip.*strmask)),0.75);
    IQR = Q3-Q1;
    clipminsm = Q1 - 1*IQR;
    clipmaxsm = Q3 + 1*IQR;
    mask_inst_clip1 = mask_inst_clip;
    mask_inst_clip1(find(smcb.*mask_inst_clip.*strmask>clipmaxsm))=0;
    mask_inst_clip1(find(smcb.*mask_inst_clip.*strmask<clipminsm))=0;
    
    cbmean1 = mean(cbmap(find(mask_inst_clip1.*strmask)));
    psmean1 = mean(psmap(find(mask_inst_clip1.*strmask)));
    cbmap = cbmap - cbmean1;
    psmap = psmap - psmean1;
    cbmean = cbmean + cbmean1;
    psmean = psmean + psmean1;
    mask_inst_clip = mask_inst_clip1;
    
    %%% get smoothed FF err %%%
    sm = fillpadsmooth(cbmap,mask_inst_clip.*strmask,50);
    N2 = sum(mask_inst_clip(:).*strmask(:));
    sprintf('%.4f,%.4f',N1/1024^2*100,N2/1024^2*100)
    %%% write the data %%%
    stackmapdat(ifield).cbmap = cbmap;
    stackmapdat(ifield).psmap = psmap;
    stackmapdat(ifield).mask_inst_clip = mask_inst_clip;
    stackmapdat(ifield).strmask = maskdat.mask(ifield).m_max(20).strmask_stack;
    stackmapdat(ifield).strnum = maskdat.mask(ifield).m_max(20).strnum_stack;
    stackmapdat(ifield).psmap_FFerr = psmap + sm;
    stackmapdat(ifield).m_min_arr = m_min_arr;
    stackmapdat(ifield).m_max_arr = m_max_arr;
    stackmapdat(ifield).cbmean = cbmean;
    stackmapdat(ifield).psmean = psmean;
    
    %%% get mask_inst_clip from strmask with limits %%%
    for m_max=17:20
        strmask = maskdat.mask(ifield).m_max(m_max).strmask_stack;
        strnum = maskdat.mask(ifield).m_max(m_max).strnum_stack;
        totmask = stackmapdat(ifield).mask_inst_clip.*strmask;
        sigmask = totmask;
        sigmask((cbmap_raw>clipmax) | (cbmap_raw<clipmin)) = 0;
        sig_sp = find((totmask-sigmask)==1);
        mask_inst_clip = stackmapdat(ifield).mask_inst_clip;
        mask_inst_clip(sig_sp)=0;

        stackmapdat(ifield).m_max(m_max).strmask = strmask;
        stackmapdat(ifield).m_max(m_max).strnum = strnum;
        stackmapdat(ifield).m_max(m_max).mask_inst_clip = mask_inst_clip;
    end
    

end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/stackmapdat',savedir),'stackmapdat');


data = zeros([5,11,1024,1024]);
for ifield=4:8
    data(ifield-3,1,:,:) = stackmapdat(ifield).rawmap;
    data(ifield-3,2,:,:) = stackmapdat(ifield).mask;
    data(ifield-3,3,:,:) = stackmapdat(ifield).DCsubmap;
    data(ifield-3,4,:,:) = stackmapdat(ifield).FF;
    data(ifield-3,5,:,:) = stackmapdat(ifield).FFunholy;
    data(ifield-3,6,:,:) = stackmapdat(ifield).map;
    data(ifield-3,7,:,:) = stackmapdat(ifield).cbmap;
    data(ifield-3,8,:,:) = stackmapdat(ifield).psmap;
    data(ifield-3,9,:,:) = stackmapdat(ifield).mask_inst_clip;
    data(ifield-3,10,:,:) = stackmapdat(ifield).strmask;
    data(ifield-3,11,:,:) = stackmapdat(ifield).strnum; 
end
save(sprintf('%s/stackmapdatarr',savedir),'data');


return


