function stack_preprocess(flight,inst,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addOptional('rmin',nan,@isnumeric);
  
  p.parse(flight,inst,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  rmin     = p.Results.rmin;
  
  clear p varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if ifield == 5
        cbmap_raw = stackmapdat(ifield).map_last10 * cal(ifield).apf2nWpm2ps;
    end
    
    %%%%% PanSTARRS srcmap %%%%%%%
    srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
        num2str(inst),'/');
    psmap_raw = fits_read(strcat(srcmapdir,dt.name,'_srcmap_ps_all.fits'));
    
    %%% masks %%%
    mask_inst = stackmapdat(ifield).mask;
    strmask = maskdat.mask(ifield).strmask_stack;
    if rmin==2
        strmask = maskdat.mask(ifield).strmask_stack_rmin2;       
    end
    
    totmask = mask_inst.*strmask;

    %%% sigma clip and mean/grad sub %%%
%     sigmask1 = sigclip_mask(cbmap_raw,totmask,3,5);
%     sigmask1 = sigclip_mask(psmap_raw,sigmask1,3,5);
%     sm = fillpadsmooth(cbmap_raw,sigmask1,2);
%     sigmask2 = sigclip_mask(sm,sigmask1,3,5);
%     sm = fillpadsmooth(psmap_raw,sigmask2,2);
%     sigmask = sigclip_mask(sm,sigmask2,3,5);
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
    
    %%% write the data %%%
    stackmapdat(ifield).cbmap = cbmap;
    stackmapdat(ifield).psmap = psmap;
    stackmapdat(ifield).mask_inst_clip = mask_inst_clip;
    if rmin==2
        stackmapdat(ifield).strmask = maskdat.mask(ifield).strmask_stack_rmin2;
        stackmapdat(ifield).strnum = maskdat.mask(ifield).strnum_stack_rmin2;
    elseif isnan(rmin)
        stackmapdat(ifield).strmask = maskdat.mask(ifield).strmask_stack;
        stackmapdat(ifield).strnum = maskdat.mask(ifield).strnum_stack;
    end
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
if rmin==2
    save(sprintf('%s/stackmapdat_rmin2',savedir),'stackmapdat');
elseif isnan(rmin)
    save(sprintf('%s/stackmapdat',savedir),'stackmapdat');
end

return


